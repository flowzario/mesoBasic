
# include "IBcapsules2D.hpp"
# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
using namespace std;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

IBcapsules2D::IBcapsules2D(const CommonParams& pin, const GetPot& input_params) : p(pin)
{

	// ---------------------------------------
	// Unpack some of the 'params' data:
	// ---------------------------------------
	
	dt = p.dt;
	np = p.np;
	rank = p.rank;		
	nbrL = p.nbrL;
	nbrR = p.nbrR;
	NX = p.NX;
	NY = p.NY;
	nx = p.nx;
	ny = p.ny;
	deli = ny + 2;
	delj = 1;
	iL = p.xOff + 1;  // lower bound (integer) of rank's domain
	iU = iL + nx - 1; // upper bound (integer) of rank's domain
	xL = double(iL);  // lower bound (double)  of rank's domain
	xU = double(iU);  // upper bound (double)  of rank's domain
	
	// ---------------------------------------
	// Set some data relating to the IBM:
	// ---------------------------------------
	
	ncaps = input_params("LBApp/ncaps",1);     // number of capsules
	nnodes = input_params("LBApp/nnodes",20);  // number of IBM nodes
	h = input_params("LBApp/h",2.0);           // half-width of each node's footprint
	kstiff = input_params("LBApp/kstiff",1.0); // stiffness parameter	
	
	xLh = xL - h;  
	xUh = xU + h; 
	wdth = int(h*2); 		

	// ---------------------------------------
	// Establish array dimensions:
	// ---------------------------------------

	for (int i=0; i<nnodes; i++) {
		x.push_back(0.0);
		y.push_back(0.0);
		x0.push_back(0.0);
		y0.push_back(0.0);
		vx.push_back(0.0);
		vy.push_back(0.0);		
		fx.push_back(0.0);
		fy.push_back(0.0);
		nn1.push_back(0);
		nn2.push_back(0);
		vxSum.push_back(0.0);
		vySum.push_back(0.0);
	}

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

IBcapsules2D::~IBcapsules2D()
{

}



// -------------------------------------------------------------------------
// Setters:
// -------------------------------------------------------------------------

void IBcapsules2D::setNodePosition(int i, double xi, double yi)
{
	x[i] = xi;
	y[i] = yi;
}

void IBcapsules2D::setNodeXPosition(int i, double xi)
{
	x[i] = xi;
}

void IBcapsules2D::setNodeYPosition(int i, double yi)
{
	y[i] = yi;
}



// -------------------------------------------------------------------------
// Getters:
// -------------------------------------------------------------------------

double IBcapsules2D::getNodeXPosition(int i)
{
	return x[i];
}

double IBcapsules2D::getNodeYPosition(int i)
{
	return y[i];
}



// -------------------------------------------------------------------------
// Initialize a single circular capsule:
// -------------------------------------------------------------------------

void IBcapsules2D::initCircularCapsule(double xo, double yo, double rad)
{
	double theta = 0.0;
	double dTheta = (2.0*3.14159265)/double(nnodes);
	for (int i=0; i<nnodes; i++) {
		x[i] = xo + rad*cos(theta);
		y[i] = yo + rad*sin(theta);
		x0[i] = x[i];
		y0[i] = y[i];
		theta += dTheta;
	}
}



// -------------------------------------------------------------------------
// Interpolate velocity from the LB grid to each of the IB nodes:
// -------------------------------------------------------------------------

/*
void IBcapsules2D::interpolateVelocity(const LBfluid2D& fl)
{
				
	// ---------------------------------------
	// loop through the IB nodes:
	// ---------------------------------------
	
	for (int n=0; n<nnodes; n++) {
		
		// -----------------------------------
		// zero node velocities:
		// -----------------------------------
		
		vx[n] = 0.0;
		vy[n] = 0.0;
		
		// -----------------------------------
		// is node footprint outside domain?
		// -----------------------------------		
		
		if (isMyNode(n) == false) continue;
		
		// -----------------------------------
		// loop over node's footprint:
		// -----------------------------------	
				
		int i0 = int(floor(x[n]));  // i-index nearest particle (rounded down)
		int j0 = int(floor(y[n]));  // j-index nearest particle (rounded down)
		
		for (int ii=0; ii<wdth; ii++) {
			int i = footprintIndex(ii,i0,wdth,NX);
			if (inMyDomain(i) == true) {
				int i_local = i - iL;
				for (int jj=0; jj<wdth; jj++) {
					int j = footprintIndex(jj,j0,wdth,NY);
					int ndx = i_local*deli + j*delj;
					double rx = x[n] - double(i);
					double ry = y[n] - double(j);
					double del = phi(rx)*phi(ry);
					vx[n] += del*fl.getU(ndx);
					vy[n] += del*fl.getV(ndx);					
				}	
			}					
		}
		
	}
		
	// ---------------------------------------
	// Sum the 'vx' and 'vy' arrays on all procs:
	// ---------------------------------------
	
	MPI::COMM_WORLD.Allreduce(&vx[0],&vxSum[0],nnodes,MPI::DOUBLE,MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&vy[0],&vySum[0],nnodes,MPI::DOUBLE,MPI::SUM);
	
	for (int i=0; i<nnodes; i++) {
		vx[i] = vxSum[i];
		vy[i] = vySum[i];
		
		if (rank==0) cout << vx[i] << " " << vy[i] << endl;
			
	}
	
	MPI::COMM_WORLD.Barrier();	
		
}
*/  


void IBcapsules2D::interpolateVelocity(const LBfluid2D& fl)
{
				
	// ---------------------------------------
	// loop through the IB nodes:
	// ---------------------------------------
	
	for (int n=0; n<nnodes; n++) {
		
		// -----------------------------------
		// zero node velocities:
		// -----------------------------------
		
		vx[n] = 0.0;
		vy[n] = 0.0;
		
		// -----------------------------------
		// get nearest point (rounded down):
		// -----------------------------------	
				
		int i0 = int(floor(x[n])) + 1;
		int j0 = int(floor(y[n])) + 1;
		
		// -----------------------------------
		// does this node overlap with me?
		// -----------------------------------
		
		if (isMyNode(i0) == false) continue;
		
		// -----------------------------------
		// loop over footprint:
		// -----------------------------------
		
		for (int i=i0; i<=i0+1; i++) {
			
			// test if this point is in my domain:
			if (inMyDomain(i) == false) continue;
			
			// get local i-index:
			int iloc = i - (iL - 1);
			
			// loop over j indices:
			for (int j=j0; j<=j0+1; j++) {
				int ndx = iloc*deli + j*delj;
				double rx = x[n] - double(i-1);
				double ry = y[n] - double(j-1);
				double del = (1.0-abs(rx))*(1.0-abs(ry));
				vx[n] += del*fl.getUStar(ndx);
				vy[n] += del*fl.getVStar(ndx);
			}
		}	
				
	}
		
	// ---------------------------------------
	// Sum the 'vx' and 'vy' arrays on all procs:
	// ---------------------------------------
	
	
	MPI::COMM_WORLD.Allreduce(&vx[0],&vxSum[0],nnodes,MPI::DOUBLE,MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&vy[0],&vySum[0],nnodes,MPI::DOUBLE,MPI::SUM);
	
	for (int i=0; i<nnodes; i++) {
		vx[i] = vxSum[i];
		vy[i] = vySum[i];
	}
	
	MPI::COMM_WORLD.Barrier();	
			
}



// -------------------------------------------------------------------------
// Extrapolate force from each IB node to the LB grid:
// -------------------------------------------------------------------------

/*
void IBcapsules2D::extrapolateForce(LBfluid2D& fl)
{
				
	// ---------------------------------------
	// loop through the IB nodes:
	// ---------------------------------------
	
	for (int n=0; n<nnodes; n++) {
				
		// -----------------------------------
		// is node footprint outside domain?
		// -----------------------------------		
		
		if (isMyNode(n) == false) continue;
		
		// -----------------------------------
		// loop over node's footprint:
		// -----------------------------------	
				
		int i0 = int(floor(x[n]));  // i-index nearest particle (rounded down)
		int j0 = int(floor(y[n]));  // j-index nearest particle (rounded down)
		
		for (int ii=0; ii<wdth; ii++) {
			int i = footprintIndex(ii,i0,wdth,NX);
			if (inMyDomain(i) == true) {
				int i_local = i - iL;
				for (int jj=0; jj<wdth; jj++) {
					int j = footprintIndex(jj,j0,wdth,NY);
					int ndx = i_local*deli + j*delj;
					double rx = x[n] - double(i);
					double ry = y[n] - double(j);
					double del = phi(abs(rx))*phi(abs(ry));
					fl.addFx(ndx,del*fx[n]);
					fl.addFy(ndx,del*fy[n]);
				}	
			}					
		}
		
	}

}
*/

void IBcapsules2D::extrapolateForce(LBfluid2D& fl)
{
				
	// ---------------------------------------
	// loop through the IB nodes:
	// ---------------------------------------
	
	for (int n=0; n<nnodes; n++) {
				
		// -----------------------------------
		// get nearest point (rounded down):
		// -----------------------------------	
				
		int i0 = int(floor(x[n])) + 1;
		int j0 = int(floor(y[n])) + 1;
		
		// -----------------------------------
		// does this node overlap with me?
		// -----------------------------------
		
		if (isMyNode(i0) == false) continue;
		
		// -----------------------------------
		// loop over footprint:
		// -----------------------------------
		
		for (int i=i0; i<=i0+1; i++) {
			
			// test if this point is in domain:
			if (inMyDomain(i) == false) continue;
			
			// get local i-index:
			int iloc = i - (iL - 1);
			
			// loop over j indices:
			for (int j=j0; j<=j0+1; j++) {
				int ndx = iloc*deli + j*delj;
				double rx = x[n] - double(i-1);
				double ry = y[n] - double(j-1);
				double del = (1.0-abs(rx))*(1.0-abs(ry));
				fl.addFx(ndx,del*fx[n]);
				fl.addFy(ndx,del*fy[n]);
			}
		}
		
	}
	
}




// -------------------------------------------------------------------------
// Compute the lattice index for an iteration through a node's footprint:
// -------------------------------------------------------------------------

int IBcapsules2D::footprintIndex(int i, int i0, int w, int N)
{
	return (i0 - w/2 + 1 + i);   // need to correct for PBC's
}



// -------------------------------------------------------------------------
// Check to see if x-position is inside rank's domain:
// -------------------------------------------------------------------------

bool IBcapsules2D::inMyDomain(int i)
{
	if (i >= iL && i <= iU) {return true;}
	else {return false;}
}



// -------------------------------------------------------------------------
// Check to see if node's footprint overlaps with rank's domain:
// -------------------------------------------------------------------------

bool IBcapsules2D::isMyNode(int i)
{	
	
	if (i >= iL-1 && i <= iU) {   // need to add PBC's
		return true;
	}
	else {
		return false;
	}		
	
	/*
	if (x[n] >= xLh && x[n] <= xUh) {
		return true;
	}
	// implement PBC for first and last domain:
	else if (rank == 0 && x[n] >= (xLh + double(NX))) {
		return true;
	}
	else if (rank == (np-1) && x[n] < h) {
		return true;
	}
	else {
		return false;
	}	
	*/
	
}



// -------------------------------------------------------------------------
// Spread function:
// -------------------------------------------------------------------------

double IBcapsules2D::phi(double r)
{
	return 0.25*(1.0 + cos(3.14159*r/2.0));
}



// -------------------------------------------------------------------------
// Update node positions based on their velocities:
// -------------------------------------------------------------------------

void IBcapsules2D::updateNodePositions()
{
	for (int i=0; i<nnodes; i++) {
		x[i] += vx[i]*dt;
		y[i] += vy[i]*dt;
	}	
}



// -------------------------------------------------------------------------
// Compute the node forces based on a membrane constitutive law:
// -------------------------------------------------------------------------

void IBcapsules2D::computeNodeForces()
{
	// use the quasi-rigid model
	for (int i=0; i<nnodes; i++) {
		fx[i] = -kstiff*(x[i] - x0[i]);
		fy[i] = -kstiff*(y[i] - y0[i]);
	}	
}



// -------------------------------------------------------------------------
// Write IBM mesh to 'vtk' file:
// -------------------------------------------------------------------------

void IBcapsules2D::writeVTKFile(std::string tagname, int tagnum)
{
	
	// -----------------------------------
	//	Define the file location and name:
	// -----------------------------------

	ofstream outfile;
	std::stringstream filenamecombine;
	filenamecombine << "vtkoutput/" << tagname << "_" << tagnum << ".vtk";
	string filename = filenamecombine.str();
	outfile.open(filename.c_str(), ios::out | ios::app);

	// -----------------------------------
	//	Write the 'vtk' file header:
	// -----------------------------------

	if (rank == 0) {
		string d = "   ";
		outfile << "# vtk DataFile Version 3.1" << endl;
		outfile << "VTK file containing IBM data" << endl;
		outfile << "ASCII" << endl;
		outfile << " " << endl;
		outfile << "DATASET POLYDATA" << endl;			
	}
	
	// -----------------------------------
	//	Write the node positions:
	// -----------------------------------

	if (rank == 0) {
		outfile << " " << endl;	
		outfile << "POINTS " << nnodes << " float" << endl;
		for (int n=0; n<nnodes; n++) {
			outfile << fixed << setprecision(3) << x[n] << "  " << y[n] << "  " << 0.0 << endl;
		}
	}

	// -----------------------------------
	//	Write lines between neighboring nodes:
	// -----------------------------------

	if (rank == 0) {
		outfile << " " << endl;	
		outfile << "LINES " << nnodes << " " << 3*nnodes << endl;
		for (int n=0; n<nnodes; n++) {
			int nplus = n+1;
			if (n == nnodes-1) nplus = 0;
			outfile << "2  " << n << "  " << nplus << endl;
		}
	}
	
	// -----------------------------------
	//	Write vertices for the nodes:
	// -----------------------------------

	if (rank == 0) {
		outfile << " " << endl;	
		outfile << "VERTICES " << nnodes << " " << 2*nnodes << endl;
		for (int n=0; n<nnodes; n++) {
			outfile << "1 " << n << endl;
		}
	}
	
	// -----------------------------------
	//	Sync the processors:
	// -----------------------------------

	MPI::COMM_WORLD.Barrier();
	
}


