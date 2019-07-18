
# ifndef TIPS4_H
# define TIPS4_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class TIPS4 : public PFBaseClass {

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

public:

    TIPS4(const CommonParams&, const GetPot&);
    ~TIPS4();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

private:

    void averageDropletSize();

};

# endif  // TIPS4_H
