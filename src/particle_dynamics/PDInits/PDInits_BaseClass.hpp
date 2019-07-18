
# ifndef PDINITSBASECLASS_H
# define PDINITSBASECLASS_H

# include <vector>
# include "../../utils/GetPot"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for initial conditions in the BD App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class PDInits_BaseClass {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of PInitCond
   // sub-classes:
   // -------------------------------------------------------------------

   static PDInits_BaseClass* PDInitFactory(const GetPot&,
                                           vector<double>&,
                                           vector<double>&,
                                           vector<double>&);

   // -------------------------------------------------------------------
   // pure virtual function:
   // -------------------------------------------------------------------

   virtual void icFunc() = 0;

};

# endif  // PDINITSBASECLASS_H
