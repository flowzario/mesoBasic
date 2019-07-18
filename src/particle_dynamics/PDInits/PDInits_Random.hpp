

# ifndef RANDOM_H
# define RANDOM_H

# include "PDInits_BaseClass.hpp"

/*
   This is an implementation of the 'PInitCond' interface class.
   It initializes a set of particles at random positions in the domain.
   It requires particles to be separated by a minimum distance; consequently,
   if you try to use a number of particles that results in a particle volume
   fraction above about 10% it will take FOREVER to initialize. So only use
   this initial condition for particle volume fractions <= 10%.
   */

class Random : public PDInits_BaseClass {

    // -------------------------------------------------------------------------
    // Private class members:
    // -------------------------------------------------------------------------

    private:

        int N;
        int nx, ny, nz;
        double dx,dy,dz;
        double Lx, Ly, Lz;
        double vscl;
        double pradii;
        string pfApp;
        bool thinFilm;
        int thickness;
        vector<double>& r;
        vector<double>& v;
        vector<double>& rad;
        double calc_separation_pbc(double, double, double);

        // -------------------------------------------------------------------------
        // Public class methods:
        // -------------------------------------------------------------------------

    public:

        Random(const GetPot& p, vector<double>&,
                vector<double>&,
                vector<double>&);
        ~Random();
        void icFunc();

};

# endif  // RANDOM_H
