# ifndef LBFLUID3D_H
# define LBFLUID3D_H

# include "../../utils/CommonParams.h"
# include "Stencil.hpp"
# include <vector>
# include <string>
# include <mpi.h>


class LBfluid3D {

private:

	const CommonParams& p;
	static int instance_count;
	int tag;
	int NX,NY,NZ;
	int nx,ny,nz;
	int gx,gy,gz;
	int nxyz,gxyz;
	int deli,delj,delk;
	int rank,np;
	int xOffset;
	int nbrL,nbrR;
	double tau;
	std::vector<double> u,v,w,r;
	std::vector<double> fx,fy,fz;
	std::vector<double> f,fstream;	

public:

	LBfluid3D(const CommonParams&);
	~LBfluid3D();
	void allocateFs(const Stencil&);
	void setTau(double);
	void setRho(int,double);
	void setFx(int,double);
	void setFy(int,double);
	void setFz(int,double);
	void setU(int,double);
	void setV(int,double);
	void setW(int,double);
	double getRho(int) const;
	double getURhoDivTau(int) const;
	double getVRhoDivTau(int) const;
	double getWRhoDivTau(int) const;
	double getRhoDivTau(int) const;
	void setFtoFeq(const Stencil&);
	void macros(const Stencil&, const bool);	
	void collideStreamUpdate(const Stencil&);
	void bounceBackWallsZdir(const Stencil&);
	void writeVTKFile(std::string,int,int,int,int);
	void ghostNodesStreaming(const Stencil&);
	void ghostNodesRho();	

private:
		
	void mpiExchange(std::vector<double>&,int,int);
	int fndx(int,int,int,int,int);

};

# endif  // LBFLUID3D_H
