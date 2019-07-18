# ifndef IBCAPSULES2D_H
# define IBCAPSULES2D_H

# include "../../utils/CommonParams.h"
# include "../../utils/GetPot"
# include "LBfluid2D.hpp"
# include <vector>
# include <string>
# include <mpi.h>


class IBcapsules2D {

private:

    const CommonParams& p;
	int rank,np;
	int nbrL,nbrR;
	int ncaps;
	int nnodes;
	int iL,iU; 
	int wdth;
	int nx,ny;
	int NX,NY;
	int deli,delj;
	double dt;
	double xL,xU,h,xLh,xUh;
	double kstiff;
	std::vector<int> nn1,nn2;
	std::vector<double> x,x0,vx,fx;
	std::vector<double> y,y0,vy,fy;
	std::vector<double> vxSum,vySum;
	
public:

	IBcapsules2D(const CommonParams&, const GetPot&);
	~IBcapsules2D();	
	void setNodePosition(int,double,double);
	void setNodeXPosition(int,double);
	void setNodeYPosition(int,double);
	double getNodeXPosition(int);
	double getNodeYPosition(int);
	void initCircularCapsule(double,double,double);
	void interpolateVelocity(const LBfluid2D&);
	void extrapolateForce(LBfluid2D&);	
	void updateNodePositions();
	void computeNodeForces();
    void writeVTKFile(std::string,int);

private:
	
	int footprintIndex(int,int,int,int);
	bool isMyNode(int);
	bool inMyDomain(int);	
	double phi(double);

};

# endif  // IBCAPSULES2D_H
