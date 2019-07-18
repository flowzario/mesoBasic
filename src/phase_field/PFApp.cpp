
# include "PFApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

PFApp::PFApp(const GetPot& input_params)
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

    ptrdiff_t locsize, locnx, offx;
    fftw_mpi_init();
    locsize = fftw_mpi_local_size_3d(p.NX,p.NY,p.NZ,MPI_COMM_WORLD,&locnx,&offx);
    p.nx = locnx;
    p.ny = p.NY;
    p.nz = p.NZ;
    p.xOff = offx;

    // ---------------------------------------
    // Create a PF object:
    // ---------------------------------------

    pf_object = PFBaseClass::PFFactory(p,input_params);


}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

PFApp::~PFApp()
{
    delete pf_object;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void PFApp::initSystem()
{
    pf_object->initPhaseField();
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void PFApp::stepForward(int step)
{

    // ----------------------------------------
    //	Set the time step:
    // ----------------------------------------

    current_step = step;
    pf_object->setTimeStep(current_step);

    // ----------------------------------------
    //	Update CH system:
    // ----------------------------------------

    pf_object->updatePhaseField();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void PFApp::writeOutput(int step)
{
    pf_object->setTimeStep(step);
    pf_object->outputPhaseField();
}
