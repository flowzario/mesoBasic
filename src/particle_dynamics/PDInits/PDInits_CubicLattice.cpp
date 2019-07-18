
# include "PDInits_CubicLattice.hpp"



// -------------------------------------------------------------------------
// Constructor...
// -------------------------------------------------------------------------

CubicLattice::CubicLattice(const GetPot& p, vector<double>& rin,
        vector<double>& vin, vector<double>& radin) :
    r(rin), v(vin), rad(radin)
{
    N = p("PDApp/N",1);
    Lx = p("Domain/nx",5);
    Ly = p("Domain/ny",5);
    Lz = p("Domain/nz",5);
    Lx *= p("Domain/dx",1.0);
    Ly *= p("Domain/dy",1.0);
    Lz *= p("Domain/dz",1.0);
    vscl = p("PDApp/initial_condition/vscl",0.0);
    rscl = p("PDApp/initial_condition/rscl",0.0);
}



// -------------------------------------------------------------------------
// Destructor...
// -------------------------------------------------------------------------

CubicLattice::~CubicLattice()
{
}



// -------------------------------------------------------------------------
// Function to calculate i.c.:
// -------------------------------------------------------------------------

void CubicLattice::icFunc()
{

    // initialize particles in a cubic array:
    int nx = pow((double)N+1,1.0/3.0); //number of particles must have integral cube root
    int ny = nx;
    int nz = nx;
    int part_index;
    double a0 = (double)Lx/(double)nx;
    double xpos,ypos,zpos;
    double r1,r2,r3;
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int k=0; k<nz; k++) {
                r1 = (double)rand()/RAND_MAX-0.5;
                r3 = (double)rand()/RAND_MAX-0.5;
                r2 = (double)rand()/RAND_MAX-0.5;
                part_index = i*ny*nz+j*nz+k;
                xpos = a0/2.0 + (double)i*a0;
                ypos = a0/2.0 + (double)j*a0;
                zpos = a0/2.0 + (double)k*a0;
                r[3*part_index+0] = xpos + rscl*r1*a0;
                r[3*part_index+1] = ypos + rscl*r2*a0;
                r[3*part_index+2] = zpos + rscl*r3*a0;

            }
        }
    }

    // initialize particle velocities:
    double vsum[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<N; i++) {
        for (int k=0; k<3; k++) {
            double r = (double)rand()/RAND_MAX;
            v[i*3+k] = vscl*(r - 0.5);
            vsum[k] += v[i*3+k];
        }
    }

    // zero the total momentum:
    if (vscl > 0.0)
    {
        for (int i=0; i<N; i++) {
            for (int k=0; k<3; k++) {
                v[i*3+k] -= vsum[k]/N;
            }
        }
    }

}
