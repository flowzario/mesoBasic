
# include "TIPSiso.hpp"
# include <iostream>
# include <fstream>



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

TIPSiso::TIPSiso(const CommonParams& pin,
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

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

TIPSiso::~TIPSiso()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void TIPSiso::initPhaseField()
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

    //	---------------------------------------
    // Output the initial configuration:
    //	---------------------------------------

    current_step = 0;
    outputPhaseField();

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void TIPSiso::updatePhaseField()
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

    c.updatePBC();
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

    mu.updatePBC();
    mob.updatePBC();    
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

}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void TIPSiso::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}
