
# ifndef PDPARTICLES_H
# define PDPARTICLES_H

# include "../utils/CommonParams.h"
# include "../utils/GetPot"
# include "PDInits/PDInits_BaseClass.hpp"
# include "PDForces/PDForces_BaseClass.hpp"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for particle-dynamics classes in the PD App.
// This class serves as an interface, and contains a factory method.
// It also contains basic elements of particle simulation (e.g. velocity-
// verlet, etc.)
// ---------------------------------------------------------------------

class PDParticles {

    public:

        PDParticles(const CommonParams&, const GetPot&);
        virtual ~PDParticles();
        void updateParticles();
        void outputParticles();
        void setTimeStep(int step) {current_step = step;}
        virtual void initParticles();

    protected:

        int N;
        int rank;
        int current_step;
        int ncell, nncells;
        int ncellx,ncelly,ncellz;
        int equilSteps;
        bool flag2D;
        bool outputForces;
        double density;
        double pradii;
        double bm_str, drag_coef;
        double dt;
        double dtover2;
        double rcut;
        double rcut2;
        double cellWidth;
        double cellWidthx,cellWidthy,cellWidthz;
        std::vector <double> r,v,f;
        std::vector <double> fbrn,fe,frp;
        std::vector <double> rad;
        std::vector <double> box;
        std::vector <double> mass;
        std::vector <double> kinEnergy;
        std::vector <int> head,list;
        std::vector <int> cellmap;
        PDInits_BaseClass* icObj;
        PDForces_BaseClass* fijObj;
        void velocityHalfKick();
        void updatePositions();
        void applyBoundaryConditions();
        void zeroForces();
        void pairwiseForces();
        virtual void auxiliaryForces();
        void writeVTKFile(string,int);
        void setupParticleCells();
        int cellIndex(int,int,int);
        double calcTotalKinEnergy();
        void writeKinEnergy(std::vector<int>,std::vector<double>);
        virtual void writeAllForces();
        void writeForce(int step,std::vector<double>&,string);

};

# endif  // PDPARTICLES_H
