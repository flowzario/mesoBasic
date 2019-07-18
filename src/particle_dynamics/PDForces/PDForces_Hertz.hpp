

# ifndef HERTZ_H
# define HERTZ_H

# include "PDForces_BaseClass.hpp"

/*
   This is an implementation of the 'PDForces_BaseClass' base class.
   It calculates an inter-particle force based on the Hertz contact model.
   */

class Hertz : public PDForces_BaseClass {

    // -------------------------------------------------------------------------
    // Private class members:
    // -------------------------------------------------------------------------

    private:

        double K;
        double rij[3];
        double rijUnit[3];
        double rijMag;
        double s2s;
        double fmag;
        double rcut,rcut2;
        double box[3];
        vector<double>& r;
        vector<double>& v;
        vector<double>& f;
        vector<double>& rad;

        // -------------------------------------------------------------------------
        // Public class methods:
        // -------------------------------------------------------------------------

    public:

        Hertz(const GetPot&, vector<double>&,
                vector<double>&,
                vector<double>&,
                vector<double>&);
        ~Hertz();
        void fijFunc(int,int);
        void equilOn();
        void equilOff();

};

# endif  // HERTZ_H
