
# include "PDApp.hpp"
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

PDApp::PDApp(const GetPot& input_params)
{

    //	---------------------------------------
    //	Assign variables from 'input_params':
    //	---------------------------------------

    p.NX = input_params("Domain/nx",1);
    p.NY = input_params("Domain/ny",1);
    p.NZ = input_params("Domain/nz",1);
    p.dx = input_params("Domain/dx",1.0);
    p.dy = input_params("Domain/dy",1.0);
    p.dz = input_params("Domain/dz",1.0);
    p.dt = input_params("Time/dt",1.0);
    p.iskip = input_params("Output/iskip",1);
    p.jskip = input_params("Output/jskip",1);
    p.kskip = input_params("Output/kskip",1);
    p.LX = p.NX*p.dx;
    p.LY = p.NY*p.dy;
    p.LZ = p.NZ*p.dz;

    // get number of particles to determin if greater than zero
    numberOfParticles = input_params("PDApp/N",0);

    //	---------------------------------------
    //	Get some MPI parameters:
    //	---------------------------------------

    p.np = MPI::COMM_WORLD.Get_size();   // # of processors
    p.rank = MPI::COMM_WORLD.Get_rank(); // my processor number

    //	---------------------------------------
    // Set dimensions:
    //	---------------------------------------

    ptrdiff_t locsize, locnx, offx;
    fftw_mpi_init();
    locsize = fftw_mpi_local_size_3d(p.NX,p.NY,p.NZ,MPI_COMM_WORLD,&locnx,&offx);
    p.nx = locnx;
    p.ny = p.NY;
    p.nz = p.NZ;
    p.xOff = offx;

    //	---------------------------------------
    //	Create a PD object:
    //	---------------------------------------

    pd_object = new PDParticles(p,input_params);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

PDApp::~PDApp()
{
    delete pd_object;
}



// -------------------------------------------------------------------------
// Initialize system:
// -------------------------------------------------------------------------

void PDApp::initSystem()
{
    pd_object->initParticles();
}



// -------------------------------------------------------------------------
// Take one step forward in time:
// -------------------------------------------------------------------------

void PDApp::stepForward(int step)
{

    // ----------------------------------------
    //	Set the time step:
    // ----------------------------------------

    current_step = step;
    pd_object->setTimeStep(current_step);

    // ----------------------------------------
    //	Update Particles system:
    // ----------------------------------------

    pd_object->updateParticles();

}



// -------------------------------------------------------------------------
// Write output:
// -------------------------------------------------------------------------

void PDApp::writeOutput(int step)
{
    pd_object->setTimeStep(step);
    if(numberOfParticles > 0)
        pd_object->outputParticles();
}
