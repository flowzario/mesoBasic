

# ifndef CUBICLATTICE_H
# define CUBICLATTICE_H

# include "PDInits_BaseClass.hpp"

/*
   This is an implementation of the 'PDInits_BaseClass' base class.
   It initializes a set of particles on a 3D cubic lattice. It
   is assume that the number of particles has an integral cube
   root i.e. (N)^(1/3) is an integer.
*/

class CubicLattice : public PDInits_BaseClass {

// -------------------------------------------------------------------------
// Private class members:
// -------------------------------------------------------------------------

private:

   int N;
   double Lx, Ly, Lz;
   double vscl;
   double rscl;
   vector<double>& r;
   vector<double>& v;
   vector<double>& rad;

// -------------------------------------------------------------------------
// Public class methods:
// -------------------------------------------------------------------------

public:

    CubicLattice(const GetPot&, vector<double>&,
                                vector<double>&,
                                vector<double>&);
    ~CubicLattice();
    void icFunc();

};

# endif  // CUBICLATTICE_H
