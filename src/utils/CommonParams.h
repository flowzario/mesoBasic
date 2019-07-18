
# ifndef COMMONPARAMS_H
# define COMMONPARAMS_H

struct CommonParams{

   int nx,ny,nz;        // number of nodes in the LOCAL domain
   int NX,NY,NZ;        // number of nodes in the GLOBAL domain
   int rank;            // processor rank
   int np;              // number of processors in the MPI communicator
   int nbrL;
   int nbrR;
   int xOff;
   int iskip;
   int jskip;
   int kskip;
   int nstep;
   int numOutputs;
   double dx,dy,dz;     // spacing between nodes
   double LX,LY,LZ;     // domain dimensions in dx,dy,dz units
   double dt;           // time step

};

# endif  // COMMONPARAMS_H
