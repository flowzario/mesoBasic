
# ifndef PDAPP_H
# define PDAPP_H

# include <mpi.h>
# include <vector>
# include "../../../../usr/include/fftw3-mpi.h"
# include "../base/MesoBase.hpp"
# include "../utils/CommonParams.h"
# include "PDParticles.hpp"


class PDApp : public MesoBase {

private:

	int current_step;
	CommonParams p;
	PDParticles* pd_object;
    int numberOfParticles; // for turning of output when N=0

public:

	PDApp(const GetPot&);
	~PDApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // PDAPP_H
