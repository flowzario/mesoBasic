
# include "scmp3D.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

scmp3D::scmp3D(const CommonParams& pin, const GetPot& input_params) : 
               p(pin), fA(p), s()
{

	// ---------------------------------------
	// set needed parameters:
	// ---------------------------------------

	deli = (p.ny+2)*(p.nz+2);
	delj = (p.nz+2);
	delk = 1;
	G = input_params("LBApp/G",-6.0);
	rhoAi = input_params("LBApp/rhoAi",0.5);
	rhoAi_noise = input_params("LBApp/rhoAi_noise",0.1);
	double tauA = input_params("LBApp/tauA",1.0);
	fA.setTau(tauA);
	
	// ---------------------------------------
	// establish stencil type:
	// ---------------------------------------

	string stype = input_params("LBApp/stencil","D3Q19");
	s.setStencil(stype);
	fA.allocateFs(s);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

scmp3D::~scmp3D()
{

}



// -------------------------------------------------------------------------
// Initialize lattice-boltzmann method:
// -------------------------------------------------------------------------

void scmp3D::initLatticeBoltzmann()
{

	// ---------------------------------------
	// initialize the system:
	// ---------------------------------------
	
	srand(time(NULL)*(p.rank+1));   // set the random seed
	
	for (int i=1; i<p.nx+1; i++) {
		for (int j=1; j<p.ny+1; j++) {	
			for (int k=1; k<p.nz+1; k++) {
				int ndx = i*deli + j*delj + k*delk;
				double rA = (double)rand()/RAND_MAX;
				// assign initial rho values:
				fA.setRho(ndx,rhoAi + rhoAi_noise*(rA-0.5));
				// assign initial velocity values:
				fA.setU(ndx,0.0);
				fA.setV(ndx,0.0);
				fA.setW(ndx,0.0);				
			}			
		}
	}
	
	fA.ghostNodesRho();
	calculateShanChenForces();
	fA.setFtoFeq(s);

}



// -------------------------------------------------------------------------
// Step forward in time the lattice-boltzmann method:
// -------------------------------------------------------------------------

void scmp3D::updateLatticeBoltzmann()
{
	
    // ---------------------------------------
    // Sync the processors:
    // ---------------------------------------

    MPI::COMM_WORLD.Barrier();
	
	// ---------------------------------------
	// Update macros:
	// ---------------------------------------
	
	bool exchangeGhostRho = true;
	fA.macros(s,exchangeGhostRho);
	
	// ---------------------------------------
	// Interfluid forces & common velocities:
	// (Shan-Chen scmp model)
	// ---------------------------------------
	
	calculateShanChenForces();
		
	// ---------------------------------------
	// collide and streaming steps:
	// ---------------------------------------

	fA.collideStreamUpdate(s);

}



// -------------------------------------------------------------------------
// Write output for the lattice-boltzmann method:
// -------------------------------------------------------------------------

void scmp3D::outputLatticeBoltzmann()
{
	int iskip = p.iskip;
	int jskip = p.jskip;
	int kskip = p.kskip;
	fA.writeVTKFile("rhoA",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Shan-Chen multi-component, multi-phase model:
// -------------------------------------------------------------------------

void scmp3D::calculateShanChenForces()
{
	
	for (int i=1; i<p.nx+1; i++) {
		for (int j=1; j<p.ny+1; j++) {
			for (int k=1; k<p.nz+1; k++) {
				
				int ndx = i*deli + j*delj + k*delk;			
			
				// get rho's & psi(rho)'s
				double rhoA = fA.getRho(ndx);
				double psi_rhoA = psi(rhoA);
			
				// sum forces from neighboring nodes
				double GsumxA = 0.0;
				double GsumyA = 0.0;
				double GsumzA = 0.0;				
				for (int n=0; n<s.nn; n++) {
					int inbr = i + s.exi[n];
					int jnbr = j + s.eyi[n];
					int knbr = k + s.ezi[n];
					int nbr  = inbr*deli + jnbr*delj + knbr*delk;
					double rhoA_nbr = fA.getRho(nbr);
					double psi_rhoA_nbr = psi(rhoA_nbr);
					GsumxA += s.wa[n]*s.ex[n]*psi_rhoA_nbr;
					GsumyA += s.wa[n]*s.ey[n]*psi_rhoA_nbr;
					GsumzA += s.wa[n]*s.ez[n]*psi_rhoA_nbr;					
				}
			
				// calculate interfluid forces
				fA.setFx(ndx,-G*psi_rhoA*GsumxA);
				fA.setFy(ndx,-G*psi_rhoA*GsumyA);
				fA.setFz(ndx,-G*psi_rhoA*GsumzA);				
			
			}			
		}
	}
	
}



// -------------------------------------------------------------------------
// Potential function of the fluid density:
// -------------------------------------------------------------------------

double scmp3D::psi(double rho)
{
   return (1.0 - exp(-rho));
}
