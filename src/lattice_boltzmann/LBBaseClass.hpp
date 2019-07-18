
# ifndef LBBASECLASS_H
# define LBBASECLASS_H

# include "../utils/CommonParams.h"
# include "../utils/GetPot"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for lattice-boltzmann classes in the LB App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class LBBaseClass {

public:

   // -------------------------------------------------------------------
   // Factory method that creates objects of sub-classes:
   // -------------------------------------------------------------------

   static LBBaseClass* LBFactory(const CommonParams&, const GetPot&);

   // -------------------------------------------------------------------
   // pure virtual functions:
   // -------------------------------------------------------------------

   virtual void initLatticeBoltzmann() = 0;
   virtual void updateLatticeBoltzmann() = 0;
   virtual void outputLatticeBoltzmann() = 0;
   virtual void setTimeStep(int) = 0;
   
   // -------------------------------------------------------------------
   // virtual destructor:
   // -------------------------------------------------------------------

   virtual ~LBBaseClass()
   {
   }

};

# endif  // LBBASECLASS_H
