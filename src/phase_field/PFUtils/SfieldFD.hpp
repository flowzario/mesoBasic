
# ifndef SFIELDFD_H
# define SFIELDFD_H

# include "../../utils/CommonParams.h"
# include <string>
# include <vector>
# include <mpi.h>


class SfieldFD {

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
        double dx2,dy2,dz2;
        std::vector<double> a;

    public:

        SfieldFD(const CommonParams&);
        ~SfieldFD();
        double getValue(int) const;
        void setValue(int,double);
        void addValue(int,double);
        void resetSfieldFD();
        void updatePBC();
        void updatePBCNoFluxZ();
		void updateBoundaries(bool,bool,bool);
        void mpiBorderExchange();
        double   Laplacian(int) const;
        SfieldFD Laplacian() const;
        SfieldFD Laplacian(const SfieldFD&) const;
        void writeVTKFile(std::string,int,int,int,int);
        SfieldFD& operator+=(const SfieldFD&);
        SfieldFD& operator+=(double);
        SfieldFD& operator-=(const SfieldFD&);
        SfieldFD& operator-=(double);
        SfieldFD& operator*=(const SfieldFD&);
        SfieldFD& operator*=(double);
        SfieldFD& operator/=(const SfieldFD&);
        SfieldFD& operator/=(double);
        SfieldFD& operator=(const SfieldFD&);
        SfieldFD  operator+(const SfieldFD&) const;
        SfieldFD  operator+(double) const;
        SfieldFD  operator-(const SfieldFD&) const;
        SfieldFD  operator-(double) const;
        SfieldFD  operator*(const SfieldFD&) const;
        SfieldFD  operator*(double) const;
        SfieldFD  operator/(const SfieldFD&) const;
        SfieldFD  operator/(double) const;

};

// Some non-member methods...
const SfieldFD operator+(double, const SfieldFD&);
const SfieldFD operator-(double, const SfieldFD&);
const SfieldFD operator*(double, const SfieldFD&);

# endif  // SFIELDFD_H
