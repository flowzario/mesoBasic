
# include "PDInits_Random.hpp"



// -------------------------------------------------------------------------
// Constructor...
// -------------------------------------------------------------------------

Random::Random(const GetPot& p, vector<double>& rin,
        vector<double>& vin, vector<double>& radin) :
    r(rin), v(vin), rad(radin)
{
    N = p("PDApp/N",1);
    nx = p("Domain/nx",5);
    ny = p("Domain/ny",5);
    nz = p("Domain/nz",5);
    dx = p("Domain/dx",1.0);
    dy = p("Domain/dy",1.0);
    dz = p("Domain/dz",1.0);
    vscl = p("PDApp/initial_condition/vscl",0.0);
    pradii = p("PDApp/pradii",1.0);

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
        Lz -= 2.0*(dz*(double)thickness + 1.5*pradii);
}



// -------------------------------------------------------------------------
// Destructor...
// -------------------------------------------------------------------------

Random::~Random()
{
}



// -------------------------------------------------------------------------
// Function to calculate i.c.:
// -------------------------------------------------------------------------

void Random::icFunc()
{
    // initialize particle positions at random locations in the domain
    for (int i=0; i<N; i++) {
        bool tooClose = true;
        double r1,r2,r3;
        while (tooClose) {
            tooClose = false;

            // get a random position
            r1 = (double)rand()/RAND_MAX*Lx;
            r2 = (double)rand()/RAND_MAX*Ly;
            if (thinFilm)
                r3 = (double)rand()/RAND_MAX*Lz + dz*(double)thickness + 1.5*pradii;
            else
                r3 = (double)rand()/RAND_MAX*Lz;

            // assign position
            r[i*3 + 0] = r1;
            r[i*3 + 1] = r2;
            r[i*3 + 2] = r3;

            // check to see if random position is too close to other positions
            for (int k = 0; k<i; k++) {
                double drx = calc_separation_pbc(r[3*i + 0],r[3*k + 0],Lx);
                double dry = calc_separation_pbc(r[3*i + 1],r[3*k + 1],Ly);
                double drz = calc_separation_pbc(r[3*i + 2],r[3*k + 2],Lz);
                double rij = sqrt(drx*drx+dry*dry+drz*drz);
                if (rij < 3.0*(rad[i]+rad[k])/2.0) 
                {
                    tooClose = true;
                    break;
                }
            }
        }
    }

    // initialize particle velocities:
    double vsum[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<N; i++) {
        for (int k=0; k<3; k++) {
            double rr = (double)rand()/RAND_MAX;
            v[i*3+k] = vscl*(rr - 0.5);
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



// -------------------------------------------------------------------------
// Separation based on PBCs:
// -------------------------------------------------------------------------

double Random::calc_separation_pbc(double ri, double rj, double L)
{
    double rij = abs(ri - rj);
    if (rij <= L/2.0) {
        return rij;
    } else
        return L - rij;
}
