
# include "LBBaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

# include "LBTypes/mcmp2D.hpp"
# include "LBTypes/mcmp3D.hpp"
# include "LBTypes/scmp3D.hpp"
# include "LBTypes/mcmp2DFilm.hpp"
# include "LBTypes/mcmp3DFilm.hpp"
# include "LBTypes/capsule2D.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

LBBaseClass* LBBaseClass::LBFactory(const CommonParams& p,
                                    const GetPot& input_params)
{

	// -----------------------------------
	// identify the requested object:
	// -----------------------------------

	string lb_type = input_params("LBApp/type","mcmp2D");

	// -----------------------------------
	// return the requested object:
	// -----------------------------------

	if (lb_type == "mcmp2D") return new mcmp2D(p,input_params);
	if (lb_type == "mcmp3D") return new mcmp3D(p,input_params);
	if (lb_type == "scmp3D") return new scmp3D(p,input_params);
	if (lb_type == "mcmp2DFilm") return new mcmp2DFilm(p,input_params);
	if (lb_type == "mcmp3DFilm") return new mcmp3DFilm(p,input_params);
	if (lb_type == "capsule2D") return new capsule2D(p,input_params);
	
	// if input file doesn't have a correct type return a nullptr
	return NULL;

}
