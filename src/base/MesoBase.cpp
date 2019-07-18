
# include "MesoBase.hpp"

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

// Lattice-Boltzmann classes:
//# include "../lattice_boltzmann/LBApp.hpp"

// Phase-Field classes:
# include "../phase_field/PFApp.hpp"

// Particle-Dynamics classes:
//# include "../particle_dynamics/PDApp.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the string 'specifier':
// {Note: all of the returnable objects inherent from 'MesoBase'}
// -------------------------------------------------------------------------

MesoBase* MesoBase::MesoObjectFactory(string specifier)
{

   // ------------------------------------------------
   // 'GetPot' object containing input parameters:
   // ------------------------------------------------

   GetPot InParams("input.dat");

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

//   if (specifier == "LBApp/") return new LBApp(InParams);
   if (specifier == "PFApp/") return new PFApp(InParams);
//   if (specifier == "PDApp/") return new PDApp(InParams);

}
