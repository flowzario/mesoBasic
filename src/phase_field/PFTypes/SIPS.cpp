/************************************************************
 * DISCLAIMER: 
 * SIPS class development continued with SIPSternary
 *
 * This class models ternary phase separation with a binary
 * system of polymer (c) and solvent and the non-solvent/solvent
 * exchange is modeled with a 1D solution to fick's second law.
 * 
 * Focus has been directed to SIPSternary.
 *
************************************************************/
# include "SIPS.hpp"
# include <iostream>
# include <fstream>



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

SIPS::SIPS(const CommonParams& pin,
                   const GetPot& input_params) : p(pin), c(p), s(p), n(p), chi(p)
{

    // ---------------------------------------
    // set needed parameters:
    // ---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    dt = p.dt;
    bx = input_params("PFApp/bx",0);
	 by = input_params("PFApp/by",1);
	 bz = input_params("PFApp/bz",1);
    NSdepth = input_params("PFApp/NSdepth",0); 
    deli = (nz+2)*(ny+2);
	 delj = (nz+2);
	 delk = 1;
    co = input_params("PFApp/co",0.1);
    chiPS = input_params("PFApp/chiPS",0.034);
    chiPN = input_params("PFApp/chiPN",2.5);
    chiCond = input_params("PFApp/chiCond",1.0);
    phiCutoff = input_params("PFApp/phiCutoff",0.5);
    M = input_params("PFApp/M",1.0);
    kap = input_params("PFApp/kap",1.0);
    N = input_params("PFApp/N",100.0);
    A = input_params("PFApp/A",1.0);
    Tinit = input_params("PFApp/Tinit",273.15);
    noiseStr = input_params("PFApp/noiseStr",0.1);
    D0 = input_params("PFApp/D0", 1.0);
    nu = input_params("PFApp/nu",1.0);
    gamma = input_params("PFApp/gamma", 10.0);
    Mweight = input_params("PFApp/Mweight",100.0);
    Mvolume = input_params("PFApp/Mvolume",0.1);
    
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

SIPS::~SIPS()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void SIPS::initPhaseField()
{

    //	---------------------------------------
    // initialize the concentration field:
    //	---------------------------------------
    int waterD = 0;                 // initialize counter
    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=1; i<nx+1; i++) {
    	  waterD = i + p.xOff - 1;    // NSdepth counter (water) (1)
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double val = 0.0;
                //double nVal = 0.0;
                //double sVal = 0.0;
                if (waterD >= NSdepth) {  // initializing polymer/solvent sol.
                    val = co + 0.1*(r-0.5);  
                    //nVal = 0.0; 
                    //sVal = 1.0 - val - nVal;
                    c.setValue(ndx,val);
                    //n.setValue(ndx,nVal);
                    //s.setValue(ndx,sVal);
                }           
                else {                    // initializing water bath
                    val = 0.0; 
                    //nVal = 1.0;
                    //sVal = 1.0 - val - nVal;
                    c.setValue(ndx,val);
                    //n.setValue(ndx,nVal);
                    //s.setValue(ndx,sVal);
                } 
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

void SIPS::updatePhaseField()
{

    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    c.updateBoundaries(bx,by,bz);
    // n.updateBoundaries(bx,by,bz);
    // s.updateBoundaries(bx,by,bz);
    MPI::COMM_WORLD.Barrier();
    // ----------------------------------------
    // initializing variables
    // ----------------------------------------
    double Dc = 1.0;      // polymer diffusion
    double Dn = 1.0;      // nonsolvent diffusion
    double Ds = 1.0;            // solvent diffusion
    double ddf_c = 1.0;
    double ddf_n = 1.0;
    double ddf_s = 1.0;
    double chiPS_diff = chiPS;
    double r = 1/N; // ratio of molar volumes of V_1/V_3 (V_n/V_p)
    double s_ratio = 1.0; // ratio of molar volumes of V_1/V_1 (V_n/V_s)
    SfieldFD mu_c(p);
    SfieldFD mob_c(p);
    // SfieldFD mu_n(p);
    // SfieldFD mob_n(p);
    // SfieldFD mu_s(p);
    // SfieldFD mob_s(p);

    int waterD = 0;
    for (int i=1; i<nx+1; i++) {
    	  waterD = i + p.xOff - 1;         
 		  double T = Tinit; 
 		  double kT = T/273.0;
 		  if (waterD >= NSdepth) {
            chiPS_diff = (chiPS-chiPN)*erf((i + p.xOff)/(2.0*sqrt(chiCond*double(current_step)*dt)))+chiPN;
        }
        else if (waterD < NSdepth) { 
            chiPS_diff = chiPN;
        } 
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                chi.setValue(ndx,chiPS_diff);
                double cc = c.getValue(ndx);
                double cc_fh = 0.0;  
                cc_fh = getFHcc(cc);	
                // ---------------------------------------
                // binary first derivative
                // ---------------------------------------
                double df_c = (log(cc_fh) + 1.0)/N - log(1.0-cc_fh) - 1.0 + chiPS_diff*(1.0-2.0*cc_fh);
					 // negative concentrations...
                if (cc <= 0.0) df_c = -1.5*A*sqrt(-cc);
                double lap_c = c.Laplacian(ndx);
                mu_c.setValue(ndx,df_c - kap*lap_c);
	 	          // --------------------------------------------------------
                // polymer self diffusion (Phillies) and 2nd Derivative FH
                // --------------------------------------------------------
		          double cc_phil = philliesDiffusion(cc);
                double ddf_c = secondDerFH(cc);		          
	             ddf_c *= kT;
                Dc = D0;
				    double Dp = Dc * exp (- gamma * pow(cc_phil,nu));	
                // ----------------------
                // mobility
                // ----------------------
                double Mc = Dp/ddf_c;
                double chiFraction = (chiPS_diff-chiPS)/(chiPN - chiPS);
                double mobScale = 0.95 + 0.25 * log(1.0225 - chiFraction);
                if (Mc > 1.0) Mc = 1.0;     // making mobility max = 1
   	          else if (Mc< 0.00001) Mc = 0.00001; // mobility min = 0.00001 till phicutoff
	             Mc *= mobScale;              
                /* if (phiCutoff < cc) Mc *= 0.001;
                if (Mc < 0.000001) Mc = 0.000001;
                if (Mc > D0) Mc = D0; 
                if (cc > phiCutoff) Mc *= 0.001;*/
                mob_c.setValue(ndx,Mc);
            }
        }
    }

    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------

    mu_c.updateBoundaries(bx,by,bz);
    mob_c.updateBoundaries(bx,by,bz);

    MPI::COMM_WORLD.Barrier();

    c += p.dt*(mu_c.Laplacian(mob_c)); 


 

    // ---------------------------------------
    // Add random fluctuations:
    // ---------------------------------------
	 waterD = 0;   
	 double c_check = 0.0;
	 double val = 0.0;
    for (int i=1; i<nx+1; i++) {
        //waterD = i + p.xOff - 1;    // counter for water depth
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double c_check = c.getValue(ndx);
                if (c_check >= phiCutoff) val = 0.0;
                else if (c_check < phiCutoff) val = noiseStr*(r-0.5);
                //if (waterD > NSdepth) c.addValue(ndx,p.dt*val); 
                // ¡¡¡nan issue!!! 
                // dont want noise added to water bath...
                c.addValue(ndx,p.dt*val);
            }
        }
    }
    
}

// ------------------------------------------------------
// flory huggins -- creating stability
// ------------------------------------------------------
double SIPS::getFHcc(double cc)
{
	double cc_fh = 0.0;
	if (cc < 0.0) { cc_fh = 0.001; }
	else if (cc > 1.0) { cc_fh = 0.999; }
	else { cc_fh = cc; }
	return cc_fh;
}


// -------------------------------------------------------
// phillies conversion
// -------------------------------------------------------
double SIPS::philliesDiffusion(double cc)
{
	double cc_phil = 0.0;
	if (cc >= 1.0) cc_phil = 1.0 * Mweight/Mvolume; // convert phi to g/L	
	else if (cc < 0.0) cc_phil = 0.000001 * Mweight/Mvolume; // convert phi to g/L 
	else { cc_phil = cc * Mweight/Mvolume; }// convert phi to g/L  
	return cc_phil;
}


// --------------------------------------------------------
// 2nd Derivative FH with phillies concentration
// --------------------------------------------------------
double SIPS::secondDerFH(double cc_fh)
{
	double ddf_c = 0.0;
	ddf_c = 0.5* (1.0/(N*cc_fh) + 1.0/(1.0-cc_fh)); 
	return ddf_c;
}


// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void SIPS::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
    // chi.writeVTKFile("chi",current_step,iskip,jskip,kskip);
    // s.writeVTKFile("s",current_step,iskip,jskip,kskip);
}
