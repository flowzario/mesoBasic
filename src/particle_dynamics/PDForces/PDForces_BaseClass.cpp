
# include "PDForces_BaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "PDForces_Hertz.hpp"
# include "PDForces_Dipole.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PDForces_BaseClass* PDForces_BaseClass::PDForcesFactory(const GetPot& p,
        vector<double>& r,
        vector<double>& v,
        vector<double>& f,
        vector<double>& rad)
{

    // -----------------------------------
    // return the requested object:
    // -----------------------------------

    string fij_type = p("PDApp/inter_particle_forces/type","Hertz");

    if (fij_type == "Hertz") return new Hertz(p,r,v,f,rad);
    if (fij_type == "Dipole") return new Dipole(p,r,v,f,rad);
    // if input file doesn't have a correct type return a nullptr
    return NULL;

}
