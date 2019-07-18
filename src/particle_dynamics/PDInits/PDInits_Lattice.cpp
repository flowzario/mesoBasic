
# include "PDInits_Lattice.hpp"



// -------------------------------------------------------------------------
// Constructor...
// -------------------------------------------------------------------------

Lattice::Lattice(const GetPot& p, vector<double>& rin,
        vector<double>& vin, vector<double>& radin) :
    r(rin), v(vin), rad(radin)
{
    N = p("PDApp/N",1);
    pradii = p("PDApp/pradii",1.0);
    nx = p("Domain/nx",5);
    ny = p("Domain/ny",5);
    nz = p("Domain/nz",5);
    dx = p("Domain/dx",1.0);
    dy = p("Domain/dy",1.0);
    dz = p("Domain/dz",1.0);
    vscl = p("PDApp/initial_condition/vscl",0.0);
    rscl = p("PDApp/initial_condition/rscl",0.0);
    Nx = p("PDApp/initial_condition/Nx",1);
    Ny = p("PDApp/initial_condition/Ny",1);
    Nz = p("PDApp/initial_condition/Nz",1);

    Lx = nx*dx;
    Ly = ny*dy;
    Lz = nz*dz;

    // check to see if this is a thin film simulation
    pfApp = p("PFApp/type","CHBD");
    thickness = p("PFApp/thickness",2);
    if (pfApp == "CHBDThinFilm")
        thinFilm = true;
    else
        thinFilm = false;
    if (thinFilm)
        Lz -= 2.0*(dz*(double)thickness + pradii);

    // check to make sure basic usage assumption is satisfied
    if (Nx*Ny*Nz != N)
    {
        throw runtime_error("For the lattice initial condition Nx*Ny*Nz must equal N!");
    }
}



// -------------------------------------------------------------------------
// Destructor...
// -------------------------------------------------------------------------

Lattice::~Lattice()
{
}



// -------------------------------------------------------------------------
// Function to calculate i.c.:
// -------------------------------------------------------------------------

void Lattice::icFunc()
{

    // initialize particles on lattice with random noise
    int part_index;
    double ax = (double)Lx/(double)Nx;
    double ay = (double)Ly/(double)Ny;
    double az = (double)Lz/(double)Nz;
    double xpos,ypos,zpos;
    double r1,r2,r3;
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            for (int k=0; k<Nz; k++) {
                r1 = (double)rand()/RAND_MAX-0.5;
                r3 = (double)rand()/RAND_MAX-0.5;
                r2 = (double)rand()/RAND_MAX-0.5;
                part_index = i*Ny*Nz+j*Nz+k;
                xpos = ax/2.0 + (double)i*ax;
                ypos = ay/2.0 + (double)j*ay;
                if (thinFilm)
                {
                    double zOffset = (double)thickness*dz + pradii;
                    zpos = az/2.0 + (double)k*az + zOffset;
                }
                else
                    zpos = az/2.0 + (double)k*az;
                r[3*part_index+0] = xpos + rscl*r1*ax;
                r[3*part_index+1] = ypos + rscl*r2*ay;
                r[3*part_index+2] = zpos + rscl*r3*az;

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
