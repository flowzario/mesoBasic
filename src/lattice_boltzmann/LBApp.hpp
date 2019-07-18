
# ifndef LBAPP_H
# define LBAPP_H

# include <mpi.h>
# include <vector>
# include "../base/MesoBase.hpp"
# include "../utils/CommonParams.h"
# include "LBBaseClass.hpp"


class LBApp : public MesoBase {
	
private:

	int current_step;
	CommonParams p;
	LBBaseClass* lb_object;

public:

	LBApp(const GetPot&);
	~LBApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // LBAPP_H
