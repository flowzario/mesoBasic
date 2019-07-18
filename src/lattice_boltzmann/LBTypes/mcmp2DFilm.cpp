
# include "mcmp2DFilm.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

mcmp2DFilm::mcmp2DFilm(const CommonParams& pin, const GetPot& input_params) : 
                       p(pin), fA(p), fB(p), s()
{

	// ---------------------------------------
	// set needed parameters:
	// ---------------------------------------

	deli = p.ny + 2;
	delj = 1;
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

	string stype = input_params("LBApp/stencil","D2Q9");
	s.setStencil(stype);
	fA.allocateFs(s);
	fB.allocateFs(s);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

mcmp2DFilm::~mcmp2DFilm()
{

}



// -------------------------------------------------------------------------
// Initialize lattice-boltzmann method:
// -------------------------------------------------------------------------

void mcmp2DFilm::initLatticeBoltzmann()
{

	// ---------------------------------------
	// initialize the system:
	// ---------------------------------------
	
	srand(time(NULL)*(p.rank+1));   // set the random seed
	
	for (int i=1; i<p.nx+1; i++) {
		for (int j=1; j<p.ny+1; j++) {			
			int ndx = i*deli + j*delj;
			double rA = (double)rand()/RAND_MAX;
			double rB = (double)rand()/RAND_MAX;
			// assign initial rho values:
			fA.setRho(ndx,rhoAi + rhoAi_noise*(rA-0.5));
			fB.setRho(ndx,rhoBi + rhoBi_noise*(rB-0.5));
			// assign initial velocity values:
			fA.setU(ndx,0.0);
			fA.setV(ndx,0.0);
			fB.setU(ndx,0.0);
			fB.setV(ndx,0.0);
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

void mcmp2DFilm::updateLatticeBoltzmann()
{
	
    // ---------------------------------------
    // Sync the processors:
    // ---------------------------------------

    MPI::COMM_WORLD.Barrier();
		
	// ---------------------------------------
	// Interfluid forces & common velocities:
	// (Shan-Chen mcmp model)
	// ---------------------------------------
	
	calculateShanChenForces();
		
	// ---------------------------------------
	// Update LBM:
	// ---------------------------------------

	int xBC = 0; 
	int yBC = 1; 
	bool exchangeGhostRho = true;
	fA.updateFluid_SC(s,xBC,yBC,exchangeGhostRho);
	fB.updateFluid_SC(s,xBC,yBC,exchangeGhostRho);

}



// -------------------------------------------------------------------------
// Write output for the lattice-boltzmann method:
// -------------------------------------------------------------------------

void mcmp2DFilm::outputLatticeBoltzmann()
{
	int iskip = p.iskip;
	int jskip = p.jskip;
	fA.writeVTKFile("rhoA",current_step,iskip,jskip);
}



// -------------------------------------------------------------------------
// Shan-Chen multi-component, multi-phase model:
// 
// Assume: y=1 and y=NY rows are walls
// 
// -------------------------------------------------------------------------

void mcmp2DFilm::calculateShanChenForces()
{
	
	for (int i=1; i<p.nx+1; i++) {
		for (int j=2; j<p.ny; j++) {    // do not include j=1 and j=ny
			
			int ndx = i*deli + j*delj;
			
			// get rho's & psi(rho)'s
			double rhoA = fA.getRho(ndx);
			double rhoB = fB.getRho(ndx);
			double psi_rhoA = psi(rhoA);
			double psi_rhoB = psi(rhoB);
			
			// sum forces from neighboring nodes
			double GsumxA = 0.0;
			double GsumyA = 0.0;
			double GsumxB = 0.0;
			double GsumyB = 0.0;			
			for (int n=0; n<s.nn; n++) {
				int inbr = i + s.exi[n];
				int jnbr = j + s.eyi[n];
				int nbr  = inbr*deli + jnbr*delj;
				
				if (jnbr == 1 || jnbr == p.ny) {
					// fluid-wall interactions					
					double psi_rhoA_nbr = psi(rhoA_sol);
					double psi_rhoB_nbr = psi(rhoB_sol);
					GsumxA += s.wa[n]*s.ex[n]*psi_rhoA_nbr;
					GsumyA += s.wa[n]*s.ey[n]*psi_rhoA_nbr;
					GsumxB += s.wa[n]*s.ex[n]*psi_rhoB_nbr;
					GsumyB += s.wa[n]*s.ey[n]*psi_rhoB_nbr;
				}
				else {
					// internal fluid-fluid interactions					
					double rhoA_nbr = fA.getRho(nbr);
					double rhoB_nbr = fB.getRho(nbr);
					double psi_rhoA_nbr = psi(rhoA_nbr);
					double psi_rhoB_nbr = psi(rhoB_nbr);
					GsumxA += s.wa[n]*s.ex[n]*psi_rhoA_nbr;
					GsumyA += s.wa[n]*s.ey[n]*psi_rhoA_nbr;
					GsumxB += s.wa[n]*s.ex[n]*psi_rhoB_nbr;
					GsumyB += s.wa[n]*s.ey[n]*psi_rhoB_nbr;
				}
				
			}
			
			// calculate interfluid forces
			fA.setFx(ndx,-G*psi_rhoA*GsumxB);
			fA.setFy(ndx,-G*psi_rhoA*GsumyB);
			fB.setFx(ndx,-G*psi_rhoB*GsumxA);
			fB.setFy(ndx,-G*psi_rhoB*GsumyA);
			
			// also, set common velocities (Shan-Chen model)
			double urdivtauA = fA.getURhoDivTau(ndx);
			double vrdivtauA = fA.getVRhoDivTau(ndx);
			double urdivtauB = fB.getURhoDivTau(ndx);
			double vrdivtauB = fB.getVRhoDivTau(ndx);
			double rdivtauA = fA.getRhoDivTau(ndx);
			double rdivtauB = fB.getRhoDivTau(ndx);			
			double uP = (urdivtauA + urdivtauB)/(rdivtauA + rdivtauB);
			double vP = (vrdivtauA + vrdivtauB)/(rdivtauA + rdivtauB);
			fA.setU(ndx,uP);  // u-prime  (common velocities)
			fA.setV(ndx,vP);  // v-prime
			fB.setU(ndx,uP);  // u-prime
			fB.setV(ndx,vP);  // v-prime
			
		}
	}
	
}



// -------------------------------------------------------------------------
// Potential function of the fluid density:
// -------------------------------------------------------------------------

double mcmp2DFilm::psi(double rho)
{
   return (1.0 - exp(-rho));
}
