
# ifndef MCMP2D_H
# define MCMP2D_H

# include "../LBBaseClass.hpp"
# include "../LBUtils/LBfluid2D.hpp"
# include "../LBUtils/Stencil.hpp"


class mcmp2D : public LBBaseClass {

private:

	const CommonParams& p;
	LBfluid2D fA,fB;
	Stencil s;  
	int current_step;
	int deli,delj;
	double G;
	double rhoAi,rhoBi;
	double rhoAi_noise,rhoBi_noise;

public:

	mcmp2D(const CommonParams&, const GetPot&);
	~mcmp2D();
	void setTimeStep(int step) {current_step = step;}
	void initLatticeBoltzmann();
	void updateLatticeBoltzmann();
	void outputLatticeBoltzmann();

private:

	void calculateShanChenForces();
	double psi(double);

};

# endif  // MCMP2D_H
