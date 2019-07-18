
# include "capsule2D.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

capsule2D::capsule2D(const CommonParams& pin, const GetPot& input_params) : 
                     p(pin), fA(p), c(p,input_params), s()
{

	// ---------------------------------------
	// set needed parameters:
	// ---------------------------------------

	deli = p.ny + 2;
	delj = 1;
	rhoAi = input_params("LBApp/rhoAi",0.5);
	rhoAi_noise = input_params("LBApp/rhoAi_noise",0.1);
	double tauA = input_params("LBApp/tauA",1.0);
	fA.setTau(tauA);
	
	// ---------------------------------------
	// establish stencil type:
	// ---------------------------------------

	string stype = input_params("LBApp/stencil","D2Q9");
	s.setStencil(stype);
	fA.allocateFs(s);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

capsule2D::~capsule2D()
{

}



// -------------------------------------------------------------------------
// Initialize lattice-boltzmann method:
// -------------------------------------------------------------------------

void capsule2D::initLatticeBoltzmann()
{

	// ---------------------------------------
	// initialize the fluid system:
	// ---------------------------------------
	
	srand(time(NULL)*(p.rank+1));   // set the random seed
	
	for (int i=1; i<p.nx+1; i++) {
		for (int j=1; j<p.ny+1; j++) {   // skip the j=1 and j=ny rows			
			int ndx = i*deli + j*delj;
			double rA = (double)rand()/RAND_MAX;
			// assign initial rho values:
			fA.setRho(ndx,rhoAi + rhoAi_noise*(rA-0.5));
			// assign initial velocity values:
			fA.setU(ndx,0.0);
			fA.setV(ndx,0.0);	
			// assign fluid forces:
			fA.setFx(ndx,0.0);
			fA.setFy(ndx,0.0);		
		}
	}
		
	// ---------------------------------------
	// initialize the capsule:
	// ---------------------------------------
	
	c.initCircularCapsule(30.0,30.0,8.0);
	
	// ---------------------------------------
	// Calculate initial fluid populations:
	// ---------------------------------------
	
	fA.setFtoFeq(s);

}



// -------------------------------------------------------------------------
// Step forward in time the lattice-boltzmann method:
// -------------------------------------------------------------------------

void capsule2D::updateLatticeBoltzmann()
{
	
    // ---------------------------------------
    // Sync the processors:
    // ---------------------------------------

    MPI::COMM_WORLD.Barrier();
	
	// ---------------------------------------
	// Compute forces in capsule membrane:
	// ---------------------------------------
	
	c.computeNodeForces();	
	
	// ---------------------------------------
	// Calculate forces on fluid lattice:
	// ---------------------------------------
	
	// zero all forces to start
	fA.zeroForces();
	
	// body force
	for (int i=1; i<p.nx+1; i++) {
		for (int j=2; j<p.ny; j++) {	// skip the j=1 and j=ny rows		
			int ndx = i*deli + j*delj;
			fA.setFx(ndx,0.000008);
			fA.setFy(ndx,0.0);	
		}
	}
	
	// add forces from capsule membrane	
	c.extrapolateForce(fA);
	
	// ---------------------------------------
	// Update LBM:
	// ---------------------------------------
	
	int xBC = 0; 
	int yBC = 1; 
	bool exchangeGhostRho = false;
	fA.updateFluid_SC(s,xBC,yBC,exchangeGhostRho);
					
	// ---------------------------------------
	// Interpolate fluid velocity to capsule:
	// ---------------------------------------

	c.interpolateVelocity(fA);
	
	// ---------------------------------------
	// Update capsule position:
	// ---------------------------------------
	
	c.updateNodePositions();
	
}



// -------------------------------------------------------------------------
// Write output for the lattice-boltzmann method:
// -------------------------------------------------------------------------

void capsule2D::outputLatticeBoltzmann()
{
	int iskip = p.iskip;
	int jskip = p.jskip;
	fA.writeVTKFile("rhoA",current_step,iskip,jskip);
	c.writeVTKFile("capsule",current_step);
}

