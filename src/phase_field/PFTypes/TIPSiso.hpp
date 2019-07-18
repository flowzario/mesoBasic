
# ifndef TIPSISO_H
# define TIPSISO_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class TIPSiso : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nx,ny,nz;
    int deli,delj,delk;
    int nxyz;
    int outAnalysisInterval;
    int numAnalysisOutputs;
    SfieldFD c;
    double co;
    double M;
    double N;
    double alpha;
    double beta;
    double kap;
    double A;
    double Tstart;
    double Tend;
    double noiseStr;

public:

    TIPSiso(const CommonParams&, const GetPot&);
    ~TIPSiso();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

};

# endif  // TIPSISO_H
