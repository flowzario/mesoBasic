
# ifndef CAPSULE2D_H
# define CAPSULE2D_H

# include "../LBBaseClass.hpp"
# include "../LBUtils/LBfluid2D.hpp"
# include "../LBUtils/IBcapsules2D.hpp"
# include "../LBUtils/Stencil.hpp"


class capsule2D : public LBBaseClass {

private:

	const CommonParams& p;
	LBfluid2D fA;
	IBcapsules2D c;
	Stencil s;  
	int current_step;
	int deli,delj;
	double rhoAi;
	double rhoAi_noise;

public:

	capsule2D(const CommonParams&, const GetPot&);
	~capsule2D();
	void setTimeStep(int step) {current_step = step;}
	void initLatticeBoltzmann();
	void updateLatticeBoltzmann();
	void outputLatticeBoltzmann();

};

# endif  // CAPSULE2D_H
