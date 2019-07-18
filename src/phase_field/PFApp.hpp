
# ifndef PFAPP_H
# define PFAPP_H

# include <mpi.h>
# include <vector>
# include <fftw3-mpi.h>
# include "../base/MesoBase.hpp"
# include "../utils/CommonParams.h"
# include "PFBaseClass.hpp"


class PFApp : public MesoBase {

private:

	int current_step;
	CommonParams p;
	PFBaseClass* pf_object;

public:

	PFApp(const GetPot&);
	~PFApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // PFAPP_H
