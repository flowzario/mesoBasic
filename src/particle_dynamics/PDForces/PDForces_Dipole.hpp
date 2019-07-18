

# ifndef DIPOLE_H
# define DIPOLE_H

# include "PDForces_BaseClass.hpp"

/*
   This is an implementation of the 'PDForces_BaseClass' base class.
   It calculates an inter-particle force based on the Dipole contact model.
   */

class Dipole : public PDForces_BaseClass {

    // -------------------------------------------------------------------------
    // Private class members:
    // -------------------------------------------------------------------------

    private:

        double eps,n;
        double Ex,Ey,Ez;
        double eqEx,eqEy,eqEz;
        double edir[3];
        double Emag;
        double radComp;
        double avgRad;
        double theta;
        double thetaComp;
        double EdotRij;
        double ScaleFactor;
        double normal[3];
        double thetaHat[3];
        double thetaMag;
        double rij[3];
        double rijMag;
        double rij2;
        double pradii;
        double rcut,rcut2;
        double box[3];
        bool thinFilm;
        vector<double>& r;
        vector<double>& v;
        vector<double>& f;
        vector<double>& rad;

    public:

        // -------------------------------------------------------------------------
        // Public class methods:
        // -------------------------------------------------------------------------

    public:

        Dipole(const GetPot&, vector<double>&,
                vector<double>&,
                vector<double>&,
                vector<double>&);
        ~Dipole();
        void fijFunc(int,int);
        void cross(double (&a)[3], double(&b)[3],double (&crs)[3]);
        void equilOn();
        void equilOff();

};

# endif  // DIPOLE_H
