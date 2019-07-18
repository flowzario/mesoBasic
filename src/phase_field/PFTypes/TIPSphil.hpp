/************************************************************
 * 
 * TIPSphil class 
 *
 * This class models binary phase separation via TIPS.
 * Temperature is isotropically lowered linearly with time.
 * Simulation time determines the quench rate, 
 * i.e. longer simulations model a slow quench rate and 
 * shorter simulations a faster  quench rate.
 *
 *
************************************************************/

# ifndef TIPSPHIL_H
# define TIPSPHIL_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class TIPSphil : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int mobStep;
    int nx,ny,nz;
    int deli,delj,delk;
    int nxyz;
    double dt;
    SfieldFD c;
    double co;
    double M;
    double N;
    double alpha;
    double beta;
    double kap;
    double A;
    double Tinit;
    double Tbath;
    double Tcrystal;
    double noiseStr;
    double nu;
    double gamma;
    double D0;
    double thermCond;
    double Mweight;
    double Mvolume;
    bool bx,by,bz;
public:

    TIPSphil(const CommonParams&, const GetPot&);
    ~TIPSphil();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // TIPSPHIL_H
