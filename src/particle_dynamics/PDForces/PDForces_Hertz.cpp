
# include "PDForces_Hertz.hpp"



// -------------------------------------------------------------------------
// Constructor...
// -------------------------------------------------------------------------

Hertz::Hertz(const GetPot& p, vector<double>& rin, vector<double>& vin,
        vector<double>& fin, vector<double>& radin) :
    r(rin), v(vin), f(fin), rad(radin)
{
    K = p("PDApp/inter_particle_forces/K",0.1);
    box[0] = p("Domain/nx",5);
    box[1] = p("Domain/ny",5);
    box[2] = p("Domain/nz",5);
    box[0] *= p("Domain/dx",1.0);
    box[1] *= p("Domain/dy",1.0);
    box[2] *= p("Domain/dz",1.0);
    rcut = p("PDApp/inter_particle_forces/rcut",4.0);
    rcut2 = rcut*rcut;
}



// -------------------------------------------------------------------------
// Destructor...
// -------------------------------------------------------------------------

Hertz::~Hertz()
{
}



// -------------------------------------------------------------------------
// Function to calculate fij:
// -------------------------------------------------------------------------

void Hertz::fijFunc(int i, int j)
{
    //compute the squared particle distance:
    double dr[3];
    double rij2 = 0.0;
    for (int k = 0; k<3; k++)
    {
        dr[k] = r[i*3+k] - r[j*3+k];
        dr[k] -= round(dr[k]/box[k])*box[k];  // <-- pbc's
        rij2 += dr[k]*dr[k];
    }
    double rij = sqrt(rij2);             // center-to-center dist.
    double s2s = rij - (rad[i]+rad[j]);  // surface-to-surface dist.
    // compute inter-particle forces when they touch:
    if (s2s < 0.0) 
    {
        double fij = 2.5*K*pow(-s2s,1.5);    // Hertz contact force
        for (int k=0; k<3; k++) {
            f[i*3+k] += fij*dr[k]/rij;
            f[j*3+k] -= fij*dr[k]/rij;
        }
    }
}


// ------------------------------------------------------------
// sets the appropriate parameters for particle equilibration.
// ------------------------------------------------------------

void Hertz::equilOn()
{
    // no parameters need to change for equilibration
}


// ------------------------------------------------------------
// resets parameters after particle equilibration.
// ------------------------------------------------------------

void Hertz::equilOff()
{
    // no parameters need to be reset after equilibration
}
