
# ifndef PDFORCESBASECLASS_H
# define PDFORCESBASECLASS_H

# include <vector>
# include "../../utils/GetPot"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for inter-particle forces in the PD App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class PDForces_BaseClass {

    public:

        // -------------------------------------------------------------------
        // Define factory method that creates objects of PDForces_BaseClass
        // sub-classes:
        // -------------------------------------------------------------------

        static PDForces_BaseClass* PDForcesFactory(const GetPot&,
                vector<double>&,
                vector<double>&,
                vector<double>&,
                vector<double>&);

        // -------------------------------------------------------------------
        // pure virtual function:
        // -------------------------------------------------------------------

        virtual void fijFunc(int,int) = 0;
        virtual void equilOff() = 0;
        virtual void equilOn() = 0;

};

# endif  // PDFORCESBASECLASS_H
