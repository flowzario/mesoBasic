/************************************************************
 * 
 * TIPSbathPHIL class 
 *
 * This class models binary phase separation via TIPS
 * with anisotropic quenching initiating at a surface
 * (x = 0) held at a constant temperature (Tbath). 
 *
 *
************************************************************/

# include "TIPSbathPHIL.hpp"
# include <iostream>
# include <fstream>



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

TIPSbathPHIL::TIPSbathPHIL(const CommonParams& pin,
                   const GetPot& input_params) : p(pin), c(p), TempOut(p)
{

    // ---------------------------------------
    // set needed parameters:
    // ---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    dt = p.dt;
    mobStep = input_params("PFApp/mobStep",1000000);
    deli = (nz+2)*(ny+2);
	 delj = (nz+2);
	 delk = 1;
	 bx = input_params("PFApp/bx",0);
	 by = input_params("PFApp/by",1);
	 bz = input_params("PFApp/bz",1);
    co = input_params("PFApp/co",0.5);
    M = input_params("PFApp/M",1.0);
    kap = input_params("PFApp/kap",1.0);
    alpha = input_params("PFApp/alpha",1.0);
    beta = input_params("PFApp/beta",1.0);
    N = input_params("PFApp/N",100.0);
    A = input_params("PFApp/A",1.0);
    Tbath = input_params("PFApp/Tbath",273.0);
    Tinit = input_params("PFApp/Tinit",273.0);
    Tcrystal = input_params("PFApp/Tcrystal",273.0);
    noiseStr = input_params("PFApp/noiseStr",0.1);
    thermCond = input_params("PFApp/thermCond", 1.0);
    nu = input_params("PFApp/nu",1.0);
    gamma = input_params("PFApp/gamma", 10.0);
    D0 = input_params("PFApp/D0", 0.0000001);
    Mweight = input_params("PFApp/Mweight",100.0);
    Mvolume = input_params("PFApp/Mvolume",0.1);
    writeTemp = input_params("PFApp/writeTemp",1);
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

TIPSbathPHIL::~TIPSbathPHIL()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void TIPSbathPHIL::initPhaseField()
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
                TempOut.setValue(ndx,460.0);
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

void TIPSbathPHIL::updatePhaseField()
{

    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    c.updateBoundaries(bx,by,bz);
    TempOut.updateBoundaries(bx,by,bz);
    MPI::COMM_WORLD.Barrier();
    double D = 0.0;
    double cc_phil = 0.001;
    double ddf = 1.0;
    SfieldFD mu(p);
    SfieldFD mob(p);
    for (int i=1; i<nx+1; i++) {
        // 1D thermal diffusion
        double T = (Tinit-Tbath)*erf((i + p.xOff)/(2.0*sqrt(thermCond*double(current_step)*dt)))+Tbath;
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
        			 double kT = T/273.0;
		          double chi = alpha/T + beta;
                int ndx = i*deli + j*delj + k*delk;
                TempOut.setValue(ndx,T);
                double cc = c.getValue(ndx);  
                // chemical potential...
                double df = (log(cc) + 1.0)/N - log(1.0-cc) - 1.0 + chi*(1.0-2.0*cc);
                df *= kT;
                if (cc <= 0.0) df = -1.5*A*sqrt(-cc);
                double lapc = c.Laplacian(ndx);
                mu.setValue(ndx,df - kap*lapc);
                // --------------------------------------------------------
                // polymer self diffusion (Phillies) and 2nd Derivative FH
                // --------------------------------------------------------
		          if (cc >= 1.0) { 	 // upper limits..
		          	cc_phil = 1.0 * Mweight/Mvolume; // convert phi to g/L 
		          	ddf = 0.5 * (1.0/(N*0.999) + 1.0/(1.0-0.999));
		          }
		          else if (cc < 0.0) { // lower limits.. 
		          	double cc_minus = 0.0001;
		          	cc_phil = cc_minus * Mweight/Mvolume; // convert phi to g/L 
		          	ddf = 0.5 * (1.0/(N*cc_minus) + 1.0/(1.0-cc_minus));  
		          }
		          else { 		         // no problems..
		            cc_phil = cc * Mweight/Mvolume; // convert phi to g/L 
		          	ddf = 0.5 * (1.0/(N*cc) + 1.0/(1.0-cc)); 
		          }
		          ddf *= kT;
                D = D0 * (T/Tinit);
                if (D > D0) D = D0;
                double Dp = D * exp (- gamma * pow(cc_phil,nu));
                // ----------------------
                // mobility
                // ----------------------
                double Mc = Dp/ddf;
                if (Mc < 0.000001) Mc = 0.000001;
                if (Mc > D0) Mc = D0;
                if (T <= Tcrystal) Mc *= 0.001; 
                mob.setValue(ndx,Mc);
            }
        }
    }
    TempOut.updateBoundaries(bx,by,bz);
    // output mobility VTK
    /*if (current_step == p.nstep) {
         int iskip = p.iskip;
         int jskip = p.jskip;
         int kskip = p.kskip;
      	TempOut.writeVTKFile("temp",current_step,iskip,jskip,kskip);
    }    */
    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------

    mu.updateBoundaries(bx,by,bz);
    mob.updateBoundaries(bx,by,bz);
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
                double val = noiseStr*(r-0.5);
                c.addValue(ndx,p.dt*val);
            }
        }
    }
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void TIPSbathPHIL::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("polymer",current_step,iskip,jskip,kskip);
    TempOut.writeVTKFile("temp",current_step,iskip,jskip,kskip);
}
