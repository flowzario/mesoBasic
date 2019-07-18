
# include "LBfluid3D.hpp"
# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
using namespace std;



// -------------------------------------------------------------------------
// static variable initialization:
// -------------------------------------------------------------------------

int LBfluid3D::instance_count = 0;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

LBfluid3D::LBfluid3D(const CommonParams& pin) : p(pin)
{

	// ---------------------------------------
	// Unpack some of the 'params' data:
	// ---------------------------------------

	rank = p.rank;
	NX = p.NX;
	NY = p.NY;
	NZ = p.NZ;
	nx = p.nx;
	ny = p.ny;
	nz = p.nz;
	nxyz = nx*ny*nz;
	xOffset = p.xOff;
	nbrL = p.nbrL;
	nbrR = p.nbrR;

	// ---------------------------------------
	// Establish array dimensions:
	// ---------------------------------------

	instance_count++;
	tag = instance_count;   // instance identifier
	gx = nx + 2;            // local x-dim. + ghost nodes
	gy = ny + 2;            // local y-dim. + ghost nodes
	gz = nz + 2;            // local z-dim. + ghost nodes
	gxyz = gx*gy*gz;        // total lattice size
	deli = gy*gz;           // index offset for neighbors in x-dim.
	delj = gz;              // index offset for neighbors in y-dim.
	delk = 1;

	for (int i=0; i<gxyz; i++) {
		r.push_back(0.0);
		u.push_back(0.0);
		v.push_back(0.0);
		w.push_back(0.0);
		fx.push_back(0.0);
		fy.push_back(0.0);
		fz.push_back(0.0);
	}

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

LBfluid3D::~LBfluid3D()
{

}



// -------------------------------------------------------------------------
// Allocate f and fstream arrays (depends on stencil type):
// -------------------------------------------------------------------------

void LBfluid3D::allocateFs(const Stencil& s)
{
	int size = gxyz*s.nn;
	for (int i=0; i<size; i++) {
		f.push_back(0.0);
		fstream.push_back(0.0);
	} 
}



// -------------------------------------------------------------------------
// Setters:
// -------------------------------------------------------------------------

void LBfluid3D::setTau(double val)
{
	tau = val;
}

void LBfluid3D::setRho(int i, double val)
{
	r[i] = val;
}

void LBfluid3D::setFx(int i, double val)
{
	fx[i] = val;
}

void LBfluid3D::setFy(int i, double val)
{
	fy[i] = val;
}

void LBfluid3D::setFz(int i, double val)
{
	fz[i] = val;
}

void LBfluid3D::setU(int i, double val)
{
	u[i] = val;
}

void LBfluid3D::setV(int i, double val)
{
	v[i] = val;
}

void LBfluid3D::setW(int i, double val)
{
	w[i] = val;
}



// -------------------------------------------------------------------------
// Getters:
// -------------------------------------------------------------------------

double LBfluid3D::getRho(int i) const
{
	return r[i];
}

double LBfluid3D::getURhoDivTau(int i) const
{
	return u[i]*r[i]/tau;
}

double LBfluid3D::getVRhoDivTau(int i) const
{
	return v[i]*r[i]/tau;
}

double LBfluid3D::getWRhoDivTau(int i) const
{
	return w[i]*r[i]/tau;
}

double LBfluid3D::getRhoDivTau(int i) const
{
	return r[i]/tau;
}



// -------------------------------------------------------------------------
// Set 'f' to 'feq'... this is ONLY done during initialization:
// -------------------------------------------------------------------------

void LBfluid3D::setFtoFeq(const Stencil& s)
{	
	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				int ndx = i*deli + j*delj + k*delk;
				double ueq = u[ndx] + tau*fx[ndx]/r[ndx];
				double veq = v[ndx] + tau*fy[ndx]/r[ndx];
				double weq = w[ndx] + tau*fz[ndx]/r[ndx];
				double uvw2 = ueq*ueq + veq*veq + weq*weq;
				for (int n=0; n<s.nn; n++) {
					int ndxn = ndx*s.nn + n;
					double evel = s.ex[n]*ueq + s.ey[n]*veq + s.ez[n]*weq;
					double feq = 1.0 + 3.0*evel + 4.5*evel*evel - 1.5*uvw2;
					feq *= r[ndx]*s.wa[n];
					f[ndxn] = feq;
				}
			}			
		}
	}
}



// -------------------------------------------------------------------------
// Update macro arrays: 'u', 'v', 'rho':
// -------------------------------------------------------------------------

void LBfluid3D::macros(const Stencil& s, const bool ghostRho)
{
    
	// ---------------------------------------
	// update macros:
	// ---------------------------------------
	
	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				int ndx = i*deli + j*delj + k*delk;
				double sum  = 0.0;
				double sumx = 0.0;
				double sumy = 0.0;
				double sumz = 0.0;
				for (int n=0; n<s.nn; n++) {
					int ndxn = ndx*s.nn + n;
					sum  += f[ndxn];
					sumx += f[ndxn]*s.ex[n];
					sumy += f[ndxn]*s.ey[n];
					sumz += f[ndxn]*s.ez[n];
				}
				r[ndx] = sum;
				u[ndx] = sumx/r[ndx];
				v[ndx] = sumy/r[ndx];
				w[ndx] = sumz/r[ndx];
			}			
		}
	}
	
	// ---------------------------------------
	// update rho on ghost nodes:
	// ---------------------------------------
	
	if (ghostRho) ghostNodesRho();
	
}



// -------------------------------------------------------------------------
// Streaming step:
// -------------------------------------------------------------------------

void LBfluid3D::collideStreamUpdate(const Stencil& s)
{

	// -----------------------------------
	// collision step:
	// -----------------------------------
	
	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				int ndx = i*deli + j*delj + k*delk;
				double ueq = u[ndx] + tau*fx[ndx]/r[ndx];
				double veq = v[ndx] + tau*fy[ndx]/r[ndx];
				double weq = w[ndx] + tau*fz[ndx]/r[ndx];
				double uvw2 = ueq*ueq + veq*veq + weq*weq;
				for (int n=0; n<s.nn; n++) {
					int ndxn = ndx*s.nn + n;
					double evel = s.ex[n]*ueq + s.ey[n]*veq + s.ez[n]*weq;
					double feq = 1.0 + 3.0*evel + 4.5*evel*evel - 1.5*uvw2;
					feq *= r[ndx]*s.wa[n];
					fstream[ndxn] = f[ndxn] - (f[ndxn] - feq)/tau;
				}
			}			
		}
	}
	
	// -----------------------------------
	// exchange ghost nodes:
	// -----------------------------------

	ghostNodesStreaming(s);

	// -----------------------------------
	// streaming step:
	// -----------------------------------

	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				int ndx = i*deli + j*delj + k*delk;
				for (int n=0; n<s.nn; n++) {
					int ndxn = (ndx)*s.nn + n;
					int inbr = i - s.exi[n];
					int jnbr = j - s.eyi[n];
					int knbr = k - s.ezi[n];
					int nbrn = (inbr*deli + jnbr*delj + knbr*delk)*s.nn + n;
					f[ndxn]  = fstream[nbrn];
				}
			}
		}
	}
	
}



// -------------------------------------------------------------------------
// Bounce-back conditions for walls located at z=1 and z=NZ.
// Assumptions: 
//  1.) streaming has already been performed, so we 
//  must retroactively implement the bounce-back conditions.
//  2.) the D3Q19 stencil is implemented as defined in the class, Stencil.
// -------------------------------------------------------------------------

void LBfluid3D::bounceBackWallsZdir(const Stencil& s)
{
			
	int nn = s.nn;
	
	for (int i=1; i<nx+1; i++) {
		for (int j=1; j<ny+1; j++) {
		
			// -----------------------------------
			// z=2 nodes:
			// -----------------------------------
				
			f[ fndx(i,j,2,5 ,nn) ] = fstream[ fndx(i,j,2,6 ,nn) ];
			f[ fndx(i,j,2,9 ,nn) ] = fstream[ fndx(i,j,2,10,nn) ];
			f[ fndx(i,j,2,11,nn) ] = fstream[ fndx(i,j,2,12,nn) ];
			f[ fndx(i,j,2,16,nn) ] = fstream[ fndx(i,j,2,15,nn) ];
			f[ fndx(i,j,2,18,nn) ] = fstream[ fndx(i,j,2,17,nn) ];
						
			// -----------------------------------
			// z=NZ-1 nodes:
			// -----------------------------------
				
			f[ fndx(i,j,nz-1,6 ,nn) ] = fstream[ fndx(i,j,nz-1,5 ,nn) ];
			f[ fndx(i,j,nz-1,10,nn) ] = fstream[ fndx(i,j,nz-1,9 ,nn) ];
			f[ fndx(i,j,nz-1,12,nn) ] = fstream[ fndx(i,j,nz-1,11,nn) ];	
			f[ fndx(i,j,nz-1,15,nn) ] = fstream[ fndx(i,j,nz-1,16,nn) ];	
			f[ fndx(i,j,nz-1,17,nn) ] = fstream[ fndx(i,j,nz-1,18,nn) ];	
			
		}						
	}
	
}



// -------------------------------------------------------------------------
// Write rho values to 'vtk' file:
// -------------------------------------------------------------------------

void LBfluid3D::writeVTKFile(std::string tagname, int tagnum,
                             int iskip, int jskip, int kskip)
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
		outfile << "VTK file containing grid data" << endl;
		outfile << "ASCII" << endl;
		outfile << " " << endl;
		outfile << "DATASET STRUCTURED_POINTS" << endl;
		outfile << "DIMENSIONS" << d << NX/iskip << d << NY/jskip << d << NZ/kskip << endl;
		outfile << "ORIGIN " << d << 0 << d << 0 << d << 0 << endl;
		outfile << "SPACING" << d << 1.0*iskip << d << 1.0*jskip << d << 1.0*kskip << endl;
		outfile << " " << endl;
		outfile << "POINT_DATA " << (NX/iskip)*(NY/jskip)*(NZ/kskip) << endl;
		outfile << "SCALARS " << tagname << " float" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
	}

	MPI::COMM_WORLD.Barrier();

	// -----------------------------------
	// Write the data:
	// NOTE: x-data increases fastest,
	//       then y-data
	// -----------------------------------

	int np = MPI::COMM_WORLD.Get_size();    // # of processors
	
	for (int k=1; k<nz+1; k+=kskip) {
		for (int j=1; j<ny+1; j+=jskip) {
			for (int p=0; p<np; p++) {
				if (p == rank) {
					for (int i=1; i<nx+1; i++) {
						int ig = i + xOffset;
						if (ig == 0 || ig%iskip == 0) {
							int ndx = j*delj + i*deli + k*delk;
							outfile << fixed << setprecision(3) << r[ndx] << endl;
						}
					}
				}
				MPI::COMM_WORLD.Barrier();
			}
		}
	}	

	// -----------------------------------
	//	Close the file:
	// -----------------------------------

	outfile.close();

}



// -------------------------------------------------------------------------
// Update ghost nodes for rho:
// -------------------------------------------------------------------------

void LBfluid3D::ghostNodesRho()
{
	
	// -----------------------------------
    // x-dir  (MPI communication)
    // -----------------------------------
	
	int size = gy*gz;     // size of communication
	mpiExchange(r,size,0);

    // -----------------------------------
    // y-dir
    // -----------------------------------

    for (int i=0; i<nx+2; i++) {
		for (int k=0; k<nz+2; k++) {
			r[i*deli + 0*delj + k*delk] = r[i*deli + ny*delj + k*delk];
			r[i*deli + (ny+1)*delj + k*delk] = r[i*deli + 1*delj + k*delk];
		}		
	}
	
    // -----------------------------------
    // z-dir
    // -----------------------------------
	
    for (int i=0; i<nx+2; i++) {
		for (int j=0; j<ny+2; j++) {
			r[i*deli + j*delj + 0*delk] = r[i*deli + j*delj + nz*delk];
			r[i*deli + j*delj + (nz+1)*delk] = r[i*deli + j*delj + 1*delk];
		}		
	}
	
}



// -------------------------------------------------------------------------
// Update ghost nodes for fstream:
// -------------------------------------------------------------------------

void LBfluid3D::ghostNodesStreaming(const Stencil& s)
{
	
	// -----------------------------------
    // x-dir  (MPI communication)
    // -----------------------------------
	
	int size = gy*gz*s.nn;  // size of communication
	mpiExchange(fstream,size,1);

    // -----------------------------------
    // y-dir
    // -----------------------------------

    for (int i=0; i<nx+2; i++) {
		for (int k=0; k<nz+2; k++) {
			for (int n=0; n<s.nn; n++) {
				fstream[ fndx(i,0   ,k,n,s.nn) ] = fstream[ fndx(i,ny,k,n,s.nn) ];
				fstream[ fndx(i,ny+1,k,n,s.nn) ] = fstream[ fndx(i,1 ,k,n,s.nn) ];
			}
		}
	}
	
    // -----------------------------------
    // z-dir
    // -----------------------------------
	
    for (int i=0; i<nx+2; i++) {
		for (int j=0; j<ny+2; j++) {
			for (int n=0; n<s.nn; n++) {
				fstream[ fndx(i,j,0   ,n,s.nn) ] = fstream[ fndx(i,j,nz,n,s.nn) ];
				fstream[ fndx(i,j,nz+1,n,s.nn) ] = fstream[ fndx(i,j,1 ,n,s.nn) ];
			}
		}
	}
	
}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors:
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void LBfluid3D::mpiExchange(std::vector<double>& a, int size, int cNum)
{
	
	MPI::Status status;
	
    // -----------------------------------
    // Send to left, Recv from right
    // -----------------------------------
	
	int ondx = 1*size;           // out index
	int indx = (nx+1)*size;	     // in index
	int stamp = tag*20 + cNum*2; // stamp is a unique int for each comm.
    
	MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrL,stamp,
	                         &a[indx],size,MPI::DOUBLE,nbrR,stamp,status);
	
	
	// -----------------------------------
	// Send to right, Recv from left
    // -----------------------------------

	ondx = nx*size;          // out index 
	indx = 0;	             // in index
	stamp += 1;              // update the stamp
	
	MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrR,stamp,
	                         &a[indx],size,MPI::DOUBLE,nbrL,stamp,status);
	
}



// -------------------------------------------------------------------------
// Index calculator for the 'f' and 'fstream' arrays:
// -------------------------------------------------------------------------

int LBfluid3D::fndx(int i, int j, int k, int n, int nn) 
{
	return (i*deli + j*delj + k*delk)*nn + n;
}
