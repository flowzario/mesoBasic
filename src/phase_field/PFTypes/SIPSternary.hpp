/************************************************************
 * DISCLAIMER:
 * SIPSternary class still under development...
 *
 * This class models ternary phase separation for a 
 * polymer (c), solvent (s), and non-solvent (n) system
 * using the Cahn-Hilliard equation and Flory-Huggins
 * thermodynamics of mixing
 *
************************************************************/
# ifndef SIPSternary_H
# define SIPSternary_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class SIPSternary : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nx,ny,nz;
    int NSdepth;
    int deli,delj,delk;
    int nxyz;
    int outAnalysisInterval;
    SfieldFD n;				 // nonsolvent concentration;  (1)
    SfieldFD c;             // polymer concentration (3)
    SfieldFD lap_polymer;
    SfieldFD lap_nonsolvent;
    double dt;
    double co,no;
    double chi23;
    double chi13;
    double M;
    double N;
    double alpha;
    double beta;
    double eta;
    double kap;
    double A;
    double Tbath;
    double Tinit;
    double noiseStr;
    double D0;
    double Mpp,Mnn;
    double Knn,Kpp;
    double v1,v2,v3;
    double Dg12;
    bool bx,by,bz;
public:

    SIPSternary(const CommonParams&, const GetPot&);
    ~SIPSternary();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // SIPSternary_H
