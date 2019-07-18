
# include "LBApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

LBApp::LBApp(const GetPot& input_params)
{

	// ---------------------------------------
	// Assign variables from 'input_params':
	// ---------------------------------------

	p.NX = input_params("Domain/nx",1);
	p.NY = input_params("Domain/ny",1);
	p.NZ = input_params("Domain/nz",1);
	p.dx = input_params("Domain/dx",1.0);
	p.dy = input_params("Domain/dy",1.0);
	p.dz = input_params("Domain/dz",1.0);
	p.dt = input_params("Time/dt",1.0);
	p.nstep = input_params("Time/nstep",1);
	p.iskip = input_params("Output/iskip",1);
	p.jskip = input_params("Output/jskip",1);
	p.kskip = input_params("Output/kskip",1);
	p.LX = p.NX*p.dx;
	p.LY = p.NY*p.dy;
	p.LZ = p.NZ*p.dz;

	// ---------------------------------------
	// Get some MPI parameters:
	// ---------------------------------------

	p.np = MPI::COMM_WORLD.Get_size();   // # of processors
	p.rank = MPI::COMM_WORLD.Get_rank(); // my processor number
	p.nbrL = (p.rank-1) + ((p.rank-1) < 0)*p.np;        // left proc. neighbor
	p.nbrR = (p.rank+1) - ((p.rank+1) > (p.np-1))*p.np; // right proc. neighbor

	// ---------------------------------------
	// Set dimensions:
	// ---------------------------------------
	
	int nxAve = int(p.NX/p.np);
	int nxExtra = p.NX%p.np;
	int nxLoc = nxAve;
	if ((p.rank+1) <= nxExtra) nxLoc = nxAve + 1;
	
	int xOff = 0;
	for (int i=0; i<p.rank; i++) {
		if (i+1 <= nxExtra) xOff += nxAve + 1;
		if (i+1  > nxExtra) xOff += nxAve;
	}

	p.nx = nxLoc;
	p.ny = p.NY;
	p.nz = p.NZ;
	p.xOff = xOff;

	// ---------------------------------------
	// Create a LB object:
	// ---------------------------------------

	lb_object = LBBaseClass::LBFactory(p,input_params);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

LBApp::~LBApp()
{
	delete lb_object;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void LBApp::initSystem()
{
	lb_object->initLatticeBoltzmann();
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void LBApp::stepForward(int step)
{

	// ----------------------------------------
	//	Set the time step:
	// ----------------------------------------

	current_step = step;
	lb_object->setTimeStep(current_step);

	// ----------------------------------------
	//	Update CH system:
	// ----------------------------------------

	lb_object->updateLatticeBoltzmann();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void LBApp::writeOutput(int step)
{
	lb_object->outputLatticeBoltzmann();
}
