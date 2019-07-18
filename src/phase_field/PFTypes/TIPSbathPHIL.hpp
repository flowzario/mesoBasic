/************************************************************
 * 
 * TIPSbathPHIL class 
 *
 * This class models binary phase separation via TIPS
 * with anisotropic quenching initiating at a surface
 * (x = 0) held at a constant temperature (Tbath).
 *
 *
************************************************************/

# ifndef TIPSBATHPHIL_H
# define TIPSBATHPHIL_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class TIPSbathPHIL : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int mobStep;
    int nx,ny,nz;
    int deli,delj,delk;
    int nxyz;
    SfieldFD c;
    SfieldFD TempOut;
    double dt;
    double co;
    double M;
    double N;
    double alpha;
    double beta;
    double kap;
    double A;
    double Tbath;
    double Tinit;
    double Tcrystal;
    double noiseStr;
    double thermCond;
    double nu;
    double gamma;
    double D0;
    double Mweight;
    double Mvolume;
    bool bx,by,bz;
    bool writeTemp;
public:

    TIPSbathPHIL(const CommonParams&, const GetPot&);
    ~TIPSbathPHIL();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // TIPSBATHPHIL_H
