
# include "SfieldFD.hpp"
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

int SfieldFD::instance_count = 0;



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

SfieldFD::SfieldFD(const CommonParams& pin) : p(pin)
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
    dx2 = p.dx*p.dx;
    dy2 = p.dy*p.dy;
    dz2 = p.dz*p.dz;

    // ---------------------------------------
    // Establish array dimensions:
    // ---------------------------------------

    instance_count++;
    tag = instance_count;   // array identifier
    gx = nx + 2;            // local x-dim. + ghost nodes
    gy = ny + 2;            // local y-dim. + ghost nodes
    gz = nz + 2;            // local z-dim. + ghost nodes
    gxyz = gx*gy*gz;        // total vector size
    deli = gz*gy;           // index offset for neighbors in x-dim.
    delj = gz;              // index offset for neighbors in y-dim.
    delk = 1;               // index offset for neighbors in z-dim.

    for (int i=0; i<gxyz; i++) a.push_back(0.0);

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

SfieldFD::~SfieldFD()
{

}



// -------------------------------------------------------------------------
// Setter, Getter, Adder, Resetter:
// -------------------------------------------------------------------------

void SfieldFD::setValue(int i, double val)
{
    a[i] = val;
}

void SfieldFD::addValue(int i, double val)
{
    a[i] += val;
}

double SfieldFD::getValue(int i) const
{
    return a[i];
}

void SfieldFD::resetSfieldFD()
{
    for (int i=0; i<gxyz; i++) a[i] = 0.0;
}



// -------------------------------------------------------------------------
// Periodic boundary conditions in all directions:
// -------------------------------------------------------------------------

void SfieldFD::updatePBC()
{

    // -----------------------------------
    // Set boundary conditions (x-dir.)
    // -----------------------------------

    mpiBorderExchange();
	


    // -----------------------------------
    // Set boundary conditions (y-dir.)
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        for (int k=1; k<nz+1; k++) {
            a[k + 0*delj + i*deli] = a[k + ny*delj + i*deli];
            a[k + (ny+1)*delj + i*deli] = a[k + 1*delj + i*deli];
        }
    }

    // -----------------------------------
    // Set boundary conditions (z-dir.)
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            a[0 + j*delj + i*deli] = a[nz + j*delj + i*deli];
            a[(nz+1) + j*delj + i*deli] = a[1 + j*delj + i*deli];
        }
    }

}



// -------------------------------------------------------------------------
// Periodic boundary conditions in x- & y-directions, no flux in z-dir:
// -------------------------------------------------------------------------

void SfieldFD::updatePBCNoFluxZ()
{

    // -----------------------------------
    // Set boundary conditions (x-dir.)
    // -----------------------------------

    mpiBorderExchange();

    // -----------------------------------
    // Set boundary conditions (y-dir.)
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        for (int k=1; k<nz+1; k++) {
            a[k + 0*delj + i*deli] = a[k + ny*delj + i*deli];
            a[k + (ny+1)*delj + i*deli] = a[k + 1*delj + i*deli];
        }
    }

    // -----------------------------------
    // Set boundary conditions (z-dir.)
    // -----------------------------------

    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            a[0 + j*delj + i*deli] = a[2 + j*delj + i*deli];
            a[(nz+1) + j*delj + i*deli] = a[(nz-1) + j*delj + i*deli];
        }
    }

}


// -------------------------------------------------------------------------
// Universal periodic boundary condition update: 0 = no flux, 1 = periodic
// -------------------------------------------------------------------------
void SfieldFD::updateBoundaries(bool bX, bool bY, bool bZ)
{
    // -----------------------------------
    // X-Direction Boundaries
    // -----------------------------------

	//Exchange boundaries between neighbor processors (MPI)
    mpiBorderExchange();

	if (!bX) {
		// Set no flux BC (x-dir.)
		//Override boundaries on end processors
		for (int j=1; j<ny+1; j++) {
			for (int k=1; k<nz+1; k++) {
				if (p.rank == 0){
					a[k + j*delj + 0*deli] = a[k + j*delj + 1*deli];
				}
				if (p.rank == p.np-1){
					a[k + j*delj + (nx+1)*deli] = a[k + j*delj + nx*deli];
				}
			}
		}
	}

    // -----------------------------------
    // Y-Direction Boundaries
    // -----------------------------------
	if (bY) {
		// Set periodic BC (y-dir.)
		for (int i=1; i<nx+1; i++) {
			for (int k=1; k<nz+1; k++) {
				a[k + 0*delj + i*deli] = a[k + ny*delj + i*deli];
				a[k + (ny+1)*delj + i*deli] = a[k + 1*delj + i*deli];
			}
		}
	} else {
		// Set no flux BC (y-dir.)
		for (int i=1; i<nx+1; i++) {
			for (int k=1; k<nz+1; k++) {
				a[k + 0*delj + i*deli] = a[k + 1*delj + i*deli];
				a[k + (ny+1)*delj + i*deli] = a[k + ny*delj + i*deli];
			}
		}
	}
    // -----------------------------------
    // Z-Direction Boundaries
    // -----------------------------------
	if (bZ) {
		// Set periodic BC (z-dir.)
		for (int i=1; i<nx+1; i++) {
			for (int j=1; j<ny+1; j++) {
				a[0 + j*delj + i*deli] = a[nz + j*delj + i*deli];
				a[(nz+1) + j*delj + i*deli] = a[1 + j*delj + i*deli];
			}
		}
	} else {
		// Set no flux BC (z-dir.)
		for (int i=1; i<nx+1; i++) {
			for (int j=1; j<ny+1; j++) {
				a[0 + j*delj + i*deli] = a[1 + j*delj + i*deli];
				a[(nz+1) + j*delj + i*deli] = a[(nz) + j*delj + i*deli];
			}
		}
	}
}



// -------------------------------------------------------------------------
// Exchange border data between neighboring processors:
// Note: 1D decomposition along the x-direction
// -------------------------------------------------------------------------

void SfieldFD::mpiBorderExchange()
{

    MPI::Status status;
    int size = deli;          // size of y-z face that will be sent

    // -----------------------------------
    // Send to left, Recv from right
    // -----------------------------------

    int stamp = tag*10;       // stamp is a unique int for each comm.
    int ondx = 1*deli;        // out index
    int indx = (nx+1)*deli;   // in  index

    MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrL,stamp,
            &a[indx],size,MPI::DOUBLE,nbrR,stamp,status);

    // -----------------------------------
    // Send to right, Recv from left
    // -----------------------------------

    stamp += 1;               // update the stamp
    ondx = nx*deli;           // out index
    indx = 0;                 // in  index

    MPI::COMM_WORLD.Sendrecv(&a[ondx],size,MPI::DOUBLE,nbrR,stamp,
            &a[indx],size,MPI::DOUBLE,nbrL,stamp,status);

}



// -------------------------------------------------------------------------
// Calculate "Laplacian" at one point and return a double:
// -------------------------------------------------------------------------

double SfieldFD::Laplacian(int ndx) const
{
    double lapx = (a[ndx+deli] + a[ndx-deli] - 2*a[ndx])/(dx2);
    double lapy = (a[ndx+delj] + a[ndx-delj] - 2*a[ndx])/(dy2);
    double lapz = (a[ndx+delk] + a[ndx-delk] - 2*a[ndx])/(dz2);
    return lapx + lapy + lapz;
}



// -------------------------------------------------------------------------
// Calculate "Laplacian" and return a SfieldFD:
// -------------------------------------------------------------------------

SfieldFD SfieldFD::Laplacian() const
{
    SfieldFD lapl(p);
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double lapx = (a[ndx+deli] + a[ndx-deli] - 2*a[ndx])/(dx2);
                double lapy = (a[ndx+delj] + a[ndx-delj] - 2*a[ndx])/(dy2);
                double lapz = (a[ndx+delk] + a[ndx-delk] - 2*a[ndx])/(dz2);
                lapl.setValue(ndx,lapx+lapy+lapz);
            }
        }
    }
    return lapl;
}



// -------------------------------------------------------------------------
// Calculate "Laplacian" w/ a non-uniform mobility field 'b':
// -------------------------------------------------------------------------

SfieldFD SfieldFD::Laplacian(const SfieldFD& b) const
{
    SfieldFD lapl(p);
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double bo = b.getValue(ndx);
                double bx1 = 2.0/(1.0/b.getValue(ndx-deli) + 1.0/bo);
                double bx2 = 2.0/(1.0/b.getValue(ndx+deli) + 1.0/bo);
                double by1 = 2.0/(1.0/b.getValue(ndx-delj) + 1.0/bo);
                double by2 = 2.0/(1.0/b.getValue(ndx+delj) + 1.0/bo);
                double bz1 = 2.0/(1.0/b.getValue(ndx-delk) + 1.0/bo);
                double bz2 = 2.0/(1.0/b.getValue(ndx+delk) + 1.0/bo);
                double lapx = (a[ndx-deli]*bx1 + a[ndx+deli]*bx2 - (bx1+bx2)*a[ndx])/dx2;
                double lapy = (a[ndx-delj]*by1 + a[ndx+delj]*by2 - (by1+by2)*a[ndx])/dy2;
                double lapz = (a[ndx-delk]*bz1 + a[ndx+delk]*bz2 - (bz1+bz2)*a[ndx])/dz2;
                lapl.setValue(ndx,lapx+lapy+lapz);
            }
        }
    }
    return lapl;
}



// -------------------------------------------------------------------------
// Write array values to 'vtk' file:
// -------------------------------------------------------------------------

void SfieldFD::writeVTKFile(std::string tagname, int tagnum,
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
    //	Write the data:
    // NOTE: x-data increases fastest,
    //       then y-data, then z-data
    // -----------------------------------

    int np = MPI::COMM_WORLD.Get_size();    // # of processors

    for (int k=1; k<nz+1; k+=kskip) {
        for (int j=1; j<ny+1; j+=jskip) {
            for (int r=0; r<np; r++) {
                if (r == rank) {
                    for (int i=1; i<nx+1; i++) {
                        int ig = i + xOffset;
                        if (ig == 0 || ig%iskip == 0) {
                            int ndx = k*delk + j*delj + i*deli;
                            outfile << fixed << setprecision(3) << a[ndx] << endl;
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
// Compound Assignment Operators:
// -------------------------------------------------------------------------

SfieldFD& SfieldFD::operator+=(const SfieldFD& rhs)
{
    for (int i=0; i<gxyz; i++) a[i] += rhs.a[i];
    return *this;
}

SfieldFD& SfieldFD::operator+=(double val)
{
    for (int i=0; i<gxyz; i++) a[i] += val;
    return *this;
}

SfieldFD& SfieldFD::operator-=(const SfieldFD& rhs)
{
    for (int i=0; i<gxyz; i++) a[i] -= rhs.a[i];
    return *this;
}

SfieldFD& SfieldFD::operator-=(double val)
{
    for (int i=0; i<gxyz; i++) a[i] -= val;
    return *this;
}

SfieldFD& SfieldFD::operator*=(const SfieldFD& rhs)
{
    for (int i=0; i<gxyz; i++) a[i] *= rhs.a[i];
    return *this;
}

SfieldFD& SfieldFD::operator*=(double val)
{
    for (int i=0; i<gxyz; i++) a[i] *= val;
    return *this;
}

SfieldFD& SfieldFD::operator/=(const SfieldFD& rhs)
{
    for (int i=0; i<gxyz; i++) {
        if (rhs.a[i] != 0.0) a[i] /= rhs.a[i];
        else a[i] = 0.0;
    }
    return *this;
}

SfieldFD& SfieldFD::operator/=(double val)
{
    for (int i=0; i<gxyz; i++) {
        if (val != 0.0) a[i] /= val;
        else a[i] = 0.0;
    }
    return *this;
}

SfieldFD& SfieldFD::operator=(const SfieldFD& rhs)
{
    for (int i=0; i<gxyz; i++) a[i] = rhs.a[i];
    return *this;
}



// -------------------------------------------------------------------------
// Binary Operators:
// -------------------------------------------------------------------------

SfieldFD SfieldFD::operator+(const SfieldFD& rhs) const
{
    SfieldFD result(p);
    result = *this;
    result += rhs;
    return result;
}

SfieldFD SfieldFD::operator+(double val) const
{
    SfieldFD result(p);
    result = *this;
    result += val;
    return result;
}

SfieldFD SfieldFD::operator-(const SfieldFD& rhs) const
{
    SfieldFD result(p);
    result = *this;
    result -= rhs;
    return result;
}

SfieldFD SfieldFD::operator-(double val) const
{
    SfieldFD result(p);
    result = *this;
    result -= val;
    return result;
}

SfieldFD SfieldFD::operator*(const SfieldFD& rhs) const
{
    SfieldFD result(p);
    result = *this;
    result *= rhs;
    return result;
}

SfieldFD SfieldFD::operator*(double val) const
{
    SfieldFD result(p);
    result = *this;
    result *= val;
    return result;
}

SfieldFD SfieldFD::operator/(const SfieldFD& rhs) const
{
    SfieldFD result(p);
    result = *this;
    result /= rhs;
    return result;
}

SfieldFD SfieldFD::operator/(double val) const
{
    SfieldFD result(p);
    result = *this;
    result /= val;
    return result;
}



// -------------------------------------------------------------------------
// Some non-member methods:
// -------------------------------------------------------------------------

const SfieldFD operator+(double b, const SfieldFD& a)
{
    return a+b;
}

const SfieldFD operator-(double b, const SfieldFD& a)
{
    return a*(-1.0) + b;
}

const SfieldFD operator*(double b, const SfieldFD& a)
{
    return a*b;
}
