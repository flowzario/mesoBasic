
# ifndef MESOBASE_H
# define MESOBASE_H

# include "../utils/GetPot"
# include <string>
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for the Meso project.  All application classes
// in the Meso project inherent from this class.
// ---------------------------------------------------------------------

class MesoBase {

public:

   // -------------------------------------------------------------------
   // Define factory method that creates objects of MesoBase sub-classes:
   // -------------------------------------------------------------------

   static MesoBase* MesoObjectFactory(string specifier);

   // -------------------------------------------------------------------
   // All sub-classes must define the below pure virtual functions:
   // -------------------------------------------------------------------

   virtual void initSystem() = 0;
   virtual void stepForward(int) = 0;
   virtual void writeOutput(int) = 0;

   // -------------------------------------------------------------------
   // Virtual destructor:
   // -------------------------------------------------------------------

   virtual ~MesoBase()
   {
   }
   
};

# endif  // MESOBASE_H
