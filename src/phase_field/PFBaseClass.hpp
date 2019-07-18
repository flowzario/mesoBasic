
# ifndef PFBASECLASS_H
# define PFBASECLASS_H

# include "../utils/CommonParams.h"
# include "../utils/GetPot"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for phase-field classes in the PF App.
// This class serves as an interface, and contains a factory method.
// ---------------------------------------------------------------------

class PFBaseClass {

public:

   // -------------------------------------------------------------------
   // Factory method that creates objects of sub-classes:
   // -------------------------------------------------------------------

   static PFBaseClass* PFFactory(const CommonParams&, const GetPot&);

   // -------------------------------------------------------------------
   // pure virtual functions:
   // -------------------------------------------------------------------

   virtual void initPhaseField() = 0;
   virtual void updatePhaseField() = 0;
   virtual void outputPhaseField() = 0;
   virtual void setTimeStep(int) = 0;

   // -------------------------------------------------------------------
   // virtual destructor:
   // -------------------------------------------------------------------

   virtual ~PFBaseClass()
   {
   }
   
};

# endif  // PFBASECLASS_H
