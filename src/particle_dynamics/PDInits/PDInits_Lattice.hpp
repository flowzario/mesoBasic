

# ifndef LATTICE_H
# define LATTICE_H

# include "PDInits_BaseClass.hpp"
# include <stdexcept>

/*
   This is an implementation of the 'PDInits_BaseClass' base class.
   It initializes a set of particles on a 3D lattice. It requires
   that the number of particles in each of the 3D directions be 
   specified in the input file as Nx, Ny, and Nz. It is also required
   that Nx*Ny*Nz = N where N is the total number of particles specified
   in the input file.
   */

class Lattice : public PDInits_BaseClass {

    // -------------------------------------------------------------------------
    // Private class members:
    // -------------------------------------------------------------------------

    private:

        int N;
        int Nx, Ny, Nz; // number of particles in each direction
        int nx, ny, nz;
        double dx,dy,dz;
        double Lx, Ly, Lz;
        double vscl;
        double rscl;
        double pradii;
        string pfApp;
        bool thinFilm;
        int thickness;
        vector<double>& r;
        vector<double>& v;
        vector<double>& rad;

        // -------------------------------------------------------------------------
        // Public class methods:
        // -------------------------------------------------------------------------

    public:

        Lattice(const GetPot&, vector<double>&,
                vector<double>&,
                vector<double>&);
        ~Lattice();
        void icFunc();

};

# endif  // LATTICE_H
