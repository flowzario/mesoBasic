
# include "PFBaseClass.hpp"
# include <string>

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

//# include "PFTypes/CHBasic.hpp"
//# include "PFTypes/CHBD.hpp"
//# include "PFTypes/CHBDThinFilm.hpp"
//# include "PFTypes/TIPS4.hpp"
//# include "PFTypes/BCPZone.hpp"
//# include "PFTypes/CHTernary.hpp"
# include "PFTypes/TIPSphil.hpp"
# include "PFTypes/TIPSbathPHIL.hpp"
# include "PFTypes/SIPS.hpp"
# include "PFTypes/SIPSternary.hpp"
//# include "PFTypes/OPFZoneTempFD.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the input file:
// -------------------------------------------------------------------------

PFBaseClass* PFBaseClass::PFFactory(const CommonParams& p,
                                    const GetPot& input_params)
{

    // -----------------------------------
    // identify the requested object:
    // -----------------------------------

    string pf_type = input_params("PFApp/type","CHBasic");

    // -----------------------------------
    // return the requested object:
    // -----------------------------------

//    if (pf_type == "CHBasic") return new CHBasic(p,input_params);
//    if (pf_type == "CHBD")    return new CHBD(p,input_params);
//    if (pf_type == "CHBDThinFilm")    return new CHBDThinFilm(p,input_params);
//    if (pf_type == "TIPS4")    return new TIPS4(p,input_params);
//    if (pf_type == "BCPZone") return new BCPZone(p,input_params);
//    if (pf_type == "CHTernary") return new CHTernary(p,input_params);
    if (pf_type == "TIPSphil") return new TIPSphil(p,input_params);
    if (pf_type == "TIPSbathPHIL") return new TIPSbathPHIL(p,input_params);
    if (pf_type == "SIPS") return new SIPS(p,input_params);
    if (pf_type == "SIPSternary") return new SIPSternary(p,input_params);
//	 if (pf_type == "OPFZoneTempFD") return new OPFZoneTempFD(p,input_params);

    // if input file doesn't have a correct type return a nullptr
    return NULL;
	
}
