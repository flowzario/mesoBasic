
# include "TIPS4.hpp"
# include <iostream>
# include <fstream>



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

TIPS4::TIPS4(const CommonParams& pin,
             const GetPot& input_params) : p(pin), c(p)
{

    // ---------------------------------------
    // set needed parameters:
    // ---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    deli = (nz+2)*(ny+2);
	delj = (nz+2);
	delk = 1;
    co = input_params("PFApp/co",0.5);
    M = input_params("PFApp/M",1.0);
    kap = input_params("PFApp/kap",1.0);
    alpha = input_params("PFApp/alpha",1.0);
    beta = input_params("PFApp/beta",1.0);
    N = input_params("PFApp/N",100.0);
    A = input_params("PFApp/A",1.0);
    Tstart = input_params("PFApp/Tstart",273.0);
    Tend = input_params("PFApp/Tend",273.0);
    numAnalysisOutputs = input_params("PFApp/numAnalysisOutputs",0);

    // ---------------------------------------
    // create analysis folder:
    // ---------------------------------------

    if (p.rank == 0) {
 	 	std::system("mkdir -p analysis");      // make analysis directory
 	 	std::system("exec rm -rf analysis/*"); // remove any existing files
 	}

    outAnalysisInterval = 0;
    if (numAnalysisOutputs != 0) outAnalysisInterval = p.nstep/numAnalysisOutputs;
    if (numAnalysisOutputs == 0) outAnalysisInterval = p.nstep+1;

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

TIPS4::~TIPS4()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void TIPS4::initPhaseField()
{

    //	---------------------------------------
    // initialize the concentration field:
    //	---------------------------------------

    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double val = co + 0.1*(r-0.5);
                c.setValue(ndx,val);
            }
        }
    }

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void TIPS4::updatePhaseField()
{

    // ---------------------------------------
    // calculate thermodynamics parameters
    // ---------------------------------------

    double T = Tstart - (Tstart-Tend)*(double(current_step)/double(p.nstep));
    double kT = T/273.0;
    double chi = alpha/T + beta;

    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    //c.updatePBC();
    c.updatePBCNoFluxZ();
    MPI::COMM_WORLD.Barrier();

    SfieldFD mu(p);
    SfieldFD mob(p);
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double cc = c.getValue(ndx);
                // chemical potential...
                double df = (log(cc) + 1.0)/N - log(1.0-cc) - 1.0 + chi*(1.0-2.0*cc);
                df *= kT;
                if (cc <= 0.0) df = -1.5*A*sqrt(-cc);
                double lapc = c.Laplacian(ndx);
                mu.setValue(ndx,df - kap*lapc);
                // mobility...
                double Mc = 1.0;
                if (cc > 0.1) Mc = 0.018/(pow(cc,1.75));
                mob.setValue(ndx,Mc);
            }
        }
    }

    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------

    //mu.updatePBC();
    //mob.updatePBC();
    mu.updatePBCNoFluxZ();
    mob.updatePBCNoFluxZ();
    MPI::COMM_WORLD.Barrier();

    c += p.dt*mu.Laplacian(mob);

    // ---------------------------------------
    // Add random fluctuations:
    // ---------------------------------------

    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double val = 0.1*(r-0.5);
                c.addValue(ndx,p.dt*val);
            }
        }
    }

    // ---------------------------------------
    // Calculate average droplet size:
    // ---------------------------------------

    if (current_step == 1 || current_step%outAnalysisInterval == 0) {
        averageDropletSize();
    }

}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void TIPS4::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Calculate averge droplet size:
// -------------------------------------------------------------------------

void TIPS4::averageDropletSize()
{

    // ---------------------------------------
    // Loop over grid to do some counting:
    // ---------------------------------------

    int ndropsTot = 0;
    int npoints = 0;

    for (int i=1; i<nx+1; i++) {

        // data for one column
        int ndrops = 0;
        bool indrop = false;

        // run along one column
        for (int j=1; j<ny+1; j++) {
            int ndx = i*deli + j*delj + 1*delk;
            double cc = c.getValue(ndx);
            // see if this point is inside a droplet:
            if (cc < 0.15) {
                if (!indrop) ndrops++;
                indrop = true;
                npoints++;
            } else {
                indrop = false;
            }
        }

        // keep running tally of # of drops
        ndropsTot += ndrops;

    }

    // ---------------------------------------
    // Collect information to rank=0:
    // ---------------------------------------

    int ndropsGlobal = 0;
    int npointsGlobal = 0;
    MPI::COMM_WORLD.Reduce(&ndropsTot,&ndropsGlobal,1,MPI_INT,MPI_SUM,0);
    MPI::COMM_WORLD.Reduce(&npoints,&npointsGlobal,1,MPI_INT,MPI_SUM,0);

    // ---------------------------------------
    // Output:
    // ---------------------------------------

    if (p.rank == 0) {
        // calculate average diameter
        int aveDropNum = ndropsGlobal/p.np;
        double aveDropSize = double(npointsGlobal)/double(aveDropNum);
        double aveDropDiam = 4*sqrt(aveDropSize)/3.14;
        if (ndropsGlobal == 0) aveDropDiam = 0.0;
        // write to file
        ofstream outfile;
        outfile.open("analysis/drop_size",std::ios_base::app);
        outfile << current_step*p.dt << " " << aveDropDiam << endl;
    }

}
