/************************************************************
 * DISCLAIMER:
 * SIPSternary class still under development...
 *
 * This class models ternary phase separation for a 
 * polymer (c), solvent (s), and non-solvent (n) phases
 * using the Cahn-Hilliard equation and Flory-Huggins
 * thermodynamics of mixing
 *
************************************************************/
# include "SIPSternary.hpp"
# include <iostream>
# include <fstream>



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

SIPSternary::SIPSternary(const CommonParams& pin,
                   const GetPot& input_params) : p(pin), c(p), n(p), lap_polymer(p), lap_nonsolvent(p)
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
    co = input_params("PFApp/co",0.20);
    no = input_params("PFApp/no",0.10);
    chi23 = input_params("PFApp/chiPS",0.034);
    chi13 = input_params("PFApp/chiPN",2.5);
    M = input_params("PFApp/M",1.0);
    kap = input_params("PFApp/kap",1.0);
    alpha = input_params("PFApp/alpha",0.3164);
    beta = input_params("PFApp/beta",0.4683);
    eta = input_params("PFApp/eta", 0.4993);
    N = input_params("PFApp/N",100.0);
    A = input_params("PFApp/A",1.0);
    Tinit = input_params("PFApp/Tinit",273.15);
    noiseStr = input_params("PFApp/noiseStr",0.1);
    D0 = input_params("PFApp/D0", 1.0);
    Knn = input_params("PFApp/Knn",1.0);
    Kpp = input_params("PFApp/Kpp",1.0);
    Mnn = input_params("PFApp/Mnn",1.0);
    Mpp = input_params("PFApp/Mpp",1.0);
    v1 = input_params("PFapp/v1",18.0);
    v2 = input_params("PFApp/v2",96.4);
    v3 = input_params("PFApp/v2",168.3);
    Dg12 = input_params("PFApp/Dg12",1.0);    
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

SIPSternary::~SIPSternary()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void SIPSternary::initPhaseField()
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
                double nVal = 0.0;
                double sVal = 0.0;
                if (waterD >= NSdepth) {  // initializing polymer/solvent sol.
                    val = co + 0.1*(r-0.5);  
                    nVal = no - 0.1*(r-0.5);
                    sVal = 1.0 - val - nVal;
                    c.setValue(ndx,val);
                    n.setValue(ndx,nVal);
                }           
                else {                    // initializing water bath
                    val = 0.0; 
                    sVal = 0.0;
                    nVal = 1.0;
                    c.setValue(ndx,val);
                    n.setValue(ndx,nVal);
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

void SIPSternary::updatePhaseField()
{

    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    c.updateBoundaries(bx,by,bz);
    n.updateBoundaries(bx,by,bz);
    MPI::COMM_WORLD.Barrier();
    SfieldFD mu_c(p);
    SfieldFD mob_c(p);
    SfieldFD mu_n(p);
    SfieldFD mob_n(p);
    double rRatio = v1/v3;
    double sRatio = v2/v3;
   /* for (int i=1; i<nx+1; i++) {         
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
	            int ndx = i*deli + j*delj + k*delk;    
               double lap_c = c.Laplacian(ndx);
               double lap_n = n.Laplacian(ndx);
               lap_polymer.setValue(ndx,lap_c);
               lap_nonsolvent.setValue(ndx,lap_n);         
	         }
	      }
	   }*/


    int waterD = 0;
    for (int i=1; i<nx+1; i++) {
    	  waterD = i + p.xOff - 1;         
 		  double T = Tinit; 
 		  double kT = T/273.0;
 		  int xOffset = i + p.xOff;
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                // ----------------------------------
                // getting ternary concentrations
                // ----------------------------------
                double cc = c.getValue(ndx);
					 double nn = n.getValue(ndx); 
					 double df_c = 0.0;
					 double df_n = 0.0;
					 // ----------------------------------------------
					 // begin flory huggins free energy of mixing:
					 // first derivative w/ respect to polymer
					 // ----------------------------------------------
					 
				    double g12 = 1.0;
					 if (cc <= 0.0) double df_c = -1.5*A*sqrt(-cc);
					 else {
					 	 if (nn <= 0.0) nn = 0.0;
					 	 double ss = 1.0 - cc - nn;
					 	 if (ss < 0.0) ss = 0.0;
						 double u1 = nn/(nn+ss);
						 double u2 = ss/(ss+nn);
						 double g12 = alpha + beta/(1.0 - eta*u2);
						 double df_c = 1.0/rRatio*(rRatio*log(cc) + rRatio - nn - sRatio*ss - rRatio*cc + 
										 (chi13*ss + sRatio*chi23*ss)*(nn+ss) - g12*nn*ss);
					 }
					 // ----------------------------------------------
					 // begin flory huggins free energy of mixing:
					 // first derivative w/ respect to solvent
					 // ----------------------------------------------
					 // reset variables
					 cc = c.getValue(ndx);
					 nn = n.getValue(ndx);
					 if (nn <= 0.0) double df_n = -1.5*A*sqrt(-nn);
					 else {
				 		 if (cc <= 0.0) cc = 0.0;
						 double ss = 1.0 - cc - nn;
						 if (ss <= 0.0) ss = 0.0;
					    // if (ss < 0) ss = 0.0;
					    // if (cc < 0) cc = 0.0;
						 double u1 = 1.0; //nn/(nn+ss);
						 double u2 = 1.0; //ss/(ss+nn);
						 double g12 = alpha + beta/(1.0 - eta*u2);
						 double df_n = log(nn) + 1.0 - nn - sRatio*nn - rRatio*cc + 
					               (g12*nn + chi13*cc)*(ss + cc) - sRatio*chi23*ss*cc -
					               ss*u1*u2*Dg12;
						 // double df_s = 1.0/sRatio*(sRatio*log(ss) + sRatio - nn - sRatio*ss - rRatio*cc + 
										      // (g12*nn + sRatio*chi23*cc)*(nn + cc) - chi13*nn*cc - nn*u1*u2*Dg12);
					  }
                 // laplacian of concentration
	         	 double laplap_c = c.Laplacian(ndx);
	         	 double laplap_n = n.Laplacian(ndx);
                mu_c.setValue(ndx,df_c - Kpp*laplap_c);// - kap*lap_c);
                mu_n.setValue(ndx,df_n - Knn*laplap_n);
                mob_c.setValue(ndx,Mpp); // should we comment this out?
                mob_n.setValue(ndx,Mnn); // constant mobility for now
             }
         }
     }

    /*for (int i=1; i<nx+1; i++) {         
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
	            int ndx = i*deli + j*delj + k*delk;
            	double df_c = mu_c.getValue(ndx);
            	double df_s = mu_s.getValue(ndx);
            	double lap_c = lap_polymer.Laplacian(ndx);
            	double lap_s = lap_solvent.Laplacian(ndx);
            	mu_c.setValue(ndx, df_c - Kpp*lap_c);
            	mu_s.setValue(ndx, df_s - Kss*lap_s);
            }
         }
      }*/
    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------
 
    mu_c.updateBoundaries(bx,by,bz);
    mob_c.updateBoundaries(bx,by,bz);
    mu_n.updateBoundaries(bx,by,bz);
    mob_n.updateBoundaries(bx,by,bz);

    MPI::COMM_WORLD.Barrier();
    
    c += p.dt*mu_c.Laplacian(mob_c); 
    n += p.dt*mu_n.Laplacian(mob_n);
 

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
                double val = noiseStr*(r-0.5);
                double nVal = val * (-1.0);
                //if (waterD > NSdepth) c.addValue(ndx,p.dt*val); 
                // ¡¡¡nan issue!!! 
                // dont want noise added to water bath...
                c.addValue(ndx,p.dt*val);
                n.addValue(ndx,p.dt*nVal);
            }
        }
    }
}

// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void SIPSternary::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
    n.writeVTKFile("n",current_step,iskip,jskip,kskip);
}


