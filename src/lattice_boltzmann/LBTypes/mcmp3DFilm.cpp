
# include "mcmp3DFilm.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

mcmp3DFilm::mcmp3DFilm(const CommonParams& pin, const GetPot& input_params) : 
                       p(pin), fA(p), fB(p), s()
{

	// ---------------------------------------
	// set needed parameters:
	// ---------------------------------------

	deli = (p.ny+2)*(p.nz+2);
	delj = (p.nz+2);
	delk = 1;
	G = input_params("LBApp/G",3.0);
	rhoA_sol = input_params("LBApp/rhoA_sol",0.5);
	rhoB_sol = input_params("LBApp/rhoB_sol",0.5);
	rhoAi = input_params("LBApp/rhoAi",0.5);
	rhoBi = input_params("LBApp/rhoBi",0.5);
	rhoAi_noise = input_params("LBApp/rhoAi_noise",0.1);
	rhoBi_noise = input_params("LBApp/rhoBi_noise",0.1);
	double tauA = input_params("LBApp/tauA",1.0);
	double tauB = input_params("LBApp/tauB",1.0);
	fA.setTau(tauA);
	fB.setTau(tauB);	
	
	// ---------------------------------------
	// establish stencil type:
	// ---------------------------------------

	string stype = input_params("LBApp/stencil","D3Q19");
	s.setStencil(stype);
	fA.allocateFs(s);
	fB.allocateFs(s);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

mcmp3DFilm::~mcmp3DFilm()
{

}



// -------------------------------------------------------------------------
// Initialize lattice-boltzmann method:
// -------------------------------------------------------------------------

void mcmp3DFilm::initLatticeBoltzmann()
{

	// ---------------------------------------
	// initialize the system:
	// ---------------------------------------
	
	srand(time(NULL)*(p.rank+1));   // set the random seed
	
	for (int i=1; i<p.nx+1; i++) {
		for (int j=1; j<p.ny+1; j++) {	
			for (int k=2; k<p.nz; k++) {   // neglect the k=1 and k=nz rows
				int ndx = i*deli + j*delj + k*delk;
				double rA = (double)rand()/RAND_MAX;
				double rB = (double)rand()/RAND_MAX;
				// assign initial rho values:
				fA.setRho(ndx,rhoAi + rhoAi_noise*(rA-0.5));
				fB.setRho(ndx,rhoBi + rhoBi_noise*(rB-0.5));
				// assign initial velocity values:
				fA.setU(ndx,0.0);
				fA.setV(ndx,0.0);
				fA.setW(ndx,0.0);
				fB.setU(ndx,0.0);
				fB.setV(ndx,0.0);
				fB.setW(ndx,0.0);
			}			
		}
	}
	
	fA.ghostNodesRho();
	fB.ghostNodesRho();
	calculateShanChenForces();
	fA.setFtoFeq(s);
	fB.setFtoFeq(s);

}



// -------------------------------------------------------------------------
// Step forward in time the lattice-boltzmann method:
// -------------------------------------------------------------------------

void mcmp3DFilm::updateLatticeBoltzmann()
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
	fB.macros(s,exchangeGhostRho);
	
	// ---------------------------------------
	// Interfluid forces & common velocities:
	// (Shan-Chen mcmp model)
	// ---------------------------------------
	
	calculateShanChenForces();
		
	// ---------------------------------------
	// collide and streaming steps:
	// ---------------------------------------

	fA.collideStreamUpdate(s);
	fB.collideStreamUpdate(s);
	
	// ---------------------------------------
	// bounce-back conditions at z=1 and z=nz:
	// ---------------------------------------
	
	fA.bounceBackWallsZdir(s);
	fB.bounceBackWallsZdir(s);

}



// -------------------------------------------------------------------------
// Write output for the lattice-boltzmann method:
// -------------------------------------------------------------------------

void mcmp3DFilm::outputLatticeBoltzmann()
{
	int iskip = p.iskip;
	int jskip = p.jskip;
	int kskip = p.kskip;
	fA.writeVTKFile("rhoA",current_step,iskip,jskip,kskip);
}



// -------------------------------------------------------------------------
// Shan-Chen multi-component, multi-phase model:
// 
// Assume: z=1 and z=NZ rows are walls
// 
// -------------------------------------------------------------------------

void mcmp3DFilm::calculateShanChenForces()
{
	
	for (int i=1; i<p.nx+1; i++) {
		for (int j=1; j<p.ny+1; j++) {
			for (int k=2; k<p.nz; k++) {   // neglect the k=1 and k=nz rows
				
				int ndx = i*deli + j*delj + k*delk;			
			
				// get rho's & psi(rho)'s
				double rhoA = fA.getRho(ndx);
				double rhoB = fB.getRho(ndx);
				double psi_rhoA = psi(rhoA);
				double psi_rhoB = psi(rhoB);
			
				// sum forces from neighboring nodes
				double GsumxA = 0.0;
				double GsumyA = 0.0;
				double GsumzA = 0.0;
				double GsumxB = 0.0;
				double GsumyB = 0.0;
				double GsumzB = 0.0;
				for (int n=0; n<s.nn; n++) {
					int inbr = i + s.exi[n];
					int jnbr = j + s.eyi[n];
					int knbr = k + s.ezi[n];
					int nbr  = inbr*deli + jnbr*delj + knbr*delk;
					
					if (knbr == 1 || knbr == p.nz) {
						// fluid-wall interactions
						double psi_rhoA_nbr = psi(rhoA_sol);
						double psi_rhoB_nbr = psi(rhoB_sol);
						GsumxA += s.wa[n]*s.ex[n]*psi_rhoA_nbr;
						GsumyA += s.wa[n]*s.ey[n]*psi_rhoA_nbr;
						GsumzA += s.wa[n]*s.ez[n]*psi_rhoA_nbr;
						GsumxB += s.wa[n]*s.ex[n]*psi_rhoB_nbr;
						GsumyB += s.wa[n]*s.ey[n]*psi_rhoB_nbr;
						GsumzB += s.wa[n]*s.ez[n]*psi_rhoB_nbr;
					}
					else {
						// internal fluid-fluid interactions
						double rhoA_nbr = fA.getRho(nbr);
						double rhoB_nbr = fB.getRho(nbr);
						double psi_rhoA_nbr = psi(rhoA_nbr);
						double psi_rhoB_nbr = psi(rhoB_nbr);
						GsumxA += s.wa[n]*s.ex[n]*psi_rhoA_nbr;
						GsumyA += s.wa[n]*s.ey[n]*psi_rhoA_nbr;
						GsumzA += s.wa[n]*s.ez[n]*psi_rhoA_nbr;
						GsumxB += s.wa[n]*s.ex[n]*psi_rhoB_nbr;
						GsumyB += s.wa[n]*s.ey[n]*psi_rhoB_nbr;
						GsumzB += s.wa[n]*s.ez[n]*psi_rhoB_nbr;
					}					
					
				}
			
				// calculate interfluid forces
				fA.setFx(ndx,-G*psi_rhoA*GsumxB);
				fA.setFy(ndx,-G*psi_rhoA*GsumyB);
				fA.setFz(ndx,-G*psi_rhoA*GsumzB);
				fB.setFx(ndx,-G*psi_rhoB*GsumxA);
				fB.setFy(ndx,-G*psi_rhoB*GsumyA);
				fB.setFz(ndx,-G*psi_rhoB*GsumzA);
			
				// also, set common velocities (Shan-Chen model)
				double urdivtauA = fA.getURhoDivTau(ndx);
				double vrdivtauA = fA.getVRhoDivTau(ndx);
				double wrdivtauA = fA.getWRhoDivTau(ndx);
				double urdivtauB = fB.getURhoDivTau(ndx);
				double vrdivtauB = fB.getVRhoDivTau(ndx);
				double wrdivtauB = fB.getWRhoDivTau(ndx);
				double rdivtauA = fA.getRhoDivTau(ndx);
				double rdivtauB = fB.getRhoDivTau(ndx);			
				double uP = (urdivtauA + urdivtauB)/(rdivtauA + rdivtauB);
				double vP = (vrdivtauA + vrdivtauB)/(rdivtauA + rdivtauB);
				double wP = (wrdivtauA + wrdivtauB)/(rdivtauA + rdivtauB);
				fA.setU(ndx,uP);  // u-prime  (common velocities)
				fA.setV(ndx,vP);  // v-prime
				fA.setW(ndx,wP);  // w-prime
				fB.setU(ndx,uP);  // u-prime
				fB.setV(ndx,vP);  // v-prime
				fB.setW(ndx,wP);  // w-prime
			
			}			
		}
	}
	
}



// -------------------------------------------------------------------------
// Potential function of the fluid density:
// -------------------------------------------------------------------------

double mcmp3DFilm::psi(double rho)
{
   return (1.0 - exp(-rho));
}
