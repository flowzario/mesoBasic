
# ifndef MESOEXECUTE_H
# define MESOEXECUTE_H

// ---------------------------------------------------------------------
// This is the Meso executioner class.  This class runs the simulation.
// ---------------------------------------------------------------------

# include "MesoBase.hpp"
# include <vector>
using namespace std;

class MesoExecute {

   private:

   vector<MesoBase*> mesoapps;

   public:

      MesoExecute();
      ~MesoExecute();
      void createMesoObjects();
      void executeMesoSimulation();

};

# endif  // MESOEXECUTE_H
