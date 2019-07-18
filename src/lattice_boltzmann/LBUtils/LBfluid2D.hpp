# ifndef LBFLUID2D_H
# define LBFLUID2D_H

# include "../../utils/CommonParams.h"
# include "Stencil.hpp"
# include <vector>
# include <string>
# include <mpi.h>


class LBfluid2D {

private:

    const CommonParams& p;
    static int instance_count;
    int tag;
    int NX,NY;
    int nx,ny;
    int gx,gy;
    int nxy,gxy;
    int deli,delj;
    int rank,np;
    int xOffset;
    int nbrL,nbrR;
	double tau;
	std::vector<double> u,v,r;
	std::vector<double> fx,fy;
	std::vector<double> f,fstream;	

public:

	LBfluid2D(const CommonParams&);
	~LBfluid2D();
	void allocateFs(const Stencil&);
	void setTau(double);
	void setRho(int,double);
    void setFx(int,double);
    void setFy(int,double);
    void addFx(int,double);
    void addFy(int,double);
	void setU(int,double);
	void setV(int,double);
	void zeroForces();
	double getU(int) const;
	double getV(int) const;
	double getUStar(int) const;
	double getVStar(int) const;
	double getRho(int) const;
	double getURhoDivTau(int) const;
	double getVRhoDivTau(int) const;
	double getRhoDivTau(int) const;
	void setFtoFeq(const Stencil&);
	void updateFluid_SC(const Stencil&, const int, const int, const bool);
	void bounceBackWallsXdir(const Stencil&);
	void bounceBackWallsYdir(const Stencil&);
    void writeVTKFile(std::string,int,int,int);
	void ghostNodesStreaming(const Stencil&);
	void ghostNodesRho();

private:
		
	void mpiExchange(std::vector<double>&,int,int);
	int fndx(int,int,int,int);

};

# endif  // LBFLUID2D_H
