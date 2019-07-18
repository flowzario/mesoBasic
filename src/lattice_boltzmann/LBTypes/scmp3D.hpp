
# ifndef SCMP3D_H
# define SCMP3D_H

# include "../LBBaseClass.hpp"
# include "../LBUtils/LBfluid3D.hpp"
# include "../LBUtils/Stencil.hpp"


class scmp3D : public LBBaseClass {

private:

	const CommonParams& p;
	LBfluid3D fA;
	Stencil s;  
	int current_step;
	int deli,delj,delk;
	double G;
	double rhoAi;
	double rhoAi_noise;

public:

	scmp3D(const CommonParams&, const GetPot&);
	~scmp3D();
	void setTimeStep(int step) {current_step = step;}
	void initLatticeBoltzmann();
	void updateLatticeBoltzmann();
	void outputLatticeBoltzmann();

private:

	void calculateShanChenForces();
	double psi(double);

};

# endif  // SCMP3D_H
