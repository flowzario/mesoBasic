
# ifndef MCMP3DFILM_H
# define MCMP3DFILM_H

# include "../LBBaseClass.hpp"
# include "../LBUtils/LBfluid3D.hpp"
# include "../LBUtils/Stencil.hpp"


class mcmp3DFilm : public LBBaseClass {

private:

	const CommonParams& p;
	LBfluid3D fA,fB;
	Stencil s;  
	int current_step;
	int deli,delj,delk;
	double G;	
	double rhoAi,rhoBi;
	double rhoAi_noise,rhoBi_noise;
	double rhoA_sol,rhoB_sol;

public:

	mcmp3DFilm(const CommonParams&, const GetPot&);
	~mcmp3DFilm();
	void setTimeStep(int step) {current_step = step;}
	void initLatticeBoltzmann();
	void updateLatticeBoltzmann();
	void outputLatticeBoltzmann();

private:

	void calculateShanChenForces();
	double psi(double);

};

# endif  // MCMP3DFILM_H
