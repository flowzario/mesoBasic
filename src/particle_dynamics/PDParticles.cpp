
# include "PDParticles.hpp"
# include <string>
# include <iomanip>
# include <fstream>
# include <string>
# include <sstream>
# include <stdlib.h>
# include <iostream>



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

PDParticles::PDParticles(const CommonParams& p, const GetPot& input_params)
{

    //	---------------------------------------
    //	Get parameters from 'CommonParams':
    //	---------------------------------------

    box.resize(3);
    box[0] = p.LX;
    box[1] = p.LY;
    box[2] = p.LZ;
    dt = p.dt;
    rank = p.rank;
    dtover2 = dt/2.0;
    current_step = 0;
    if (p.nz == 1) flag2D = true;
    if (p.nz  > 1) flag2D = false;

    //	---------------------------------------
    //	Get other parameters from 'GetPot':
    //	---------------------------------------

    N = input_params("PDApp/N",1);   // # of particles
    pradii = input_params("PDApp/pradii",1.0);
    rcut = input_params("PDApp/inter_particle_forces/rcut",6.0*pradii);
    rcut2 = rcut*rcut;
    density = input_params("PDApp/density",0.23873241463);
    drag_coef = input_params("PDApp/drag_coef",0.0);
    bm_str = input_params("PDApp/bm_str",0.0);
    outputForces = input_params("PDApp/outputForces",0);
    equilSteps = input_params("PDApp/equilSteps",0);

    //	---------------------------------------
    //	Establish vector dimensions:
    //	---------------------------------------

    for (int i=0; i<N; i++) 
        rad.push_back(pradii);
    double vol;
    for (int i=0; i<N; i++) {
        vol = 4.0/3.0*3.14159265359*rad[i]*rad[i]*rad[i];
        mass.push_back(density*vol);
        for (int k=0; k<3; k++) {
            r.push_back(0.0);
            v.push_back(0.0);
            f.push_back(0.0);
            fbrn.push_back(0.0);
            fe.push_back(0.0);
            frp.push_back(0.0);
        }
    }

    // ---------------------------------------
    // Create objects for init. cond. &
    // inter-particle force:
    // ---------------------------------------

    icObj = PDInits_BaseClass::PDInitFactory(input_params,r,v,rad);
    fijObj = PDForces_BaseClass::PDForcesFactory(input_params,r,v,f,rad);

    // ---------------------------------------
    // initialize linked-list cells:
    // ---------------------------------------

    setupParticleCells();

}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

PDParticles::~PDParticles()
{

}



// -------------------------------------------------------------------------
// Initializer:
// -------------------------------------------------------------------------

void PDParticles::initParticles()
{
    icObj->icFunc();

    // --------------------------------
    // equilibrate particles
    // --------------------------------

    std::vector<int> steps;
    double tke = calcTotalKinEnergy();
    kinEnergy.push_back(tke);
    steps.push_back(0);

    // turn on equilibration particle interactions
    fijObj->equilOn();

    // perform the equilibration steps
    for (int i=1; i<=equilSteps;i++)
    {
        updateParticles();
        tke = calcTotalKinEnergy();
        kinEnergy.push_back(tke);
        steps.push_back(i);
    }

    // turn off equilibration particle interactions
    fijObj->equilOff();

    // write the equilibration data to file
    writeKinEnergy(steps,kinEnergy);
    current_step = 0;
}



// -------------------------------------------------------------------------
// Updater:
// -------------------------------------------------------------------------

void PDParticles::updateParticles()
{
    zeroForces();
    pairwiseForces();
    auxiliaryForces();
    velocityHalfKick();
    velocityHalfKick();
    updatePositions();
    applyBoundaryConditions();
}



// -------------------------------------------------------------------------
// Update the particle positions by the time step:
// -------------------------------------------------------------------------

void PDParticles::updatePositions()
{
    for (int i=0; i<N; i++) {
        for (int k=0; k<3; k++) r[i*3+k] += dt*v[i*3+k];
    }
}



// -------------------------------------------------------------------------
// Update particle velocities by half the time step:
// -------------------------------------------------------------------------

void PDParticles::velocityHalfKick()
{
    for (int i=0; i<N; i++) {
        for (int k=0; k<3; k++) v[i*3+k] += dtover2*f[i*3+k]/mass[i];
    }
}



// -------------------------------------------------------------------------
// Enforce boundary conditions (periodic):
// -------------------------------------------------------------------------

void PDParticles::applyBoundaryConditions()
{
    for (int i=0; i<N; i++) {
        for (int k=0; k<3; k++) r[i*3+k] -= floor(r[i*3+k]/box[k])*box[k];
    }
}




// -------------------------------------------------------------------------
// Zero all forces:
// -------------------------------------------------------------------------

void PDParticles::zeroForces()
{
    for (int i=0; i<N; i++) {
        for (int k=0; k<3; k++) f[i*3+k] = 0.0;
    }
}




// -------------------------------------------------------------------------
// Pairwise particle forces:
// -------------------------------------------------------------------------

void PDParticles::pairwiseForces()
{

    //	---------------------------------------
    //	First, update linked-list vectors:
    //	---------------------------------------

    for (int i=0; i<ncell; i++) head[i] = -1;

    for (int i=0; i<N; i++) {
        int icell = int(floor(r[i*3+0]/cellWidthx))*ncellz*ncelly +
            int(floor(r[i*3+1]/cellWidthy))*ncellz +
            int(floor(r[i*3+2]/cellWidthz));
        list[i] = head[icell];
        head[icell] = i;
    }

    //	---------------------------------------
    //	Then, loop over cells to calculate
    // interactions:
    //	---------------------------------------

    for (int ic=0; ic<ncellx; ic++)
        for (int jc=0; jc<ncelly; jc++)
            for (int kc=0; kc<ncellz; kc++) {

                // loop over particles in current cell:
                int icell = ic*ncellz*ncelly + jc*ncellz + kc;
                int i = head[icell];
                while (i >= 0) {

                    // loop over other particles in current cell:
                    int j = list[i];
                    while (j >= 0) {
                        fijObj->fijFunc(i,j);
                        j = list[j];
                    }

                    // loop over neighboring cells:
                    for (int nbor=0; nbor<nncells; nbor++) {
                        int jcell = cellmap[icell*nncells + nbor];
                        // loop over all particles in jcell:
                        int k = head[jcell];
                        while (k >= 0) {
                            fijObj->fijFunc(i,k);
                            k = list[k];
                        }
                    }

                    // move to next particle in icell:
                    i = list[i];

                }

            }
}



// -------------------------------------------------------------------------
// Auxiliary forces:
// -------------------------------------------------------------------------

void PDParticles::auxiliaryForces()
{
    for (int i=0; i<N; i++) {
        for (int j=0; j<3; j++) {
            double rr = (double)rand()/RAND_MAX;
            f[i*3+j] += sqrt(bm_str*drag_coef*2.0)*2.0*(rr-0.5);
            f[i*3+j] -= drag_coef*v[i*3+j];
        }
    }
}



// -------------------------------------------------------------------------
// Outputer:
// -------------------------------------------------------------------------

void PDParticles::outputParticles()
{
    writeVTKFile("particles",current_step);    
    if (outputForces)
        writeAllForces();
}



// -------------------------------------------------------------------------
// Write VTK file for particles:
// -------------------------------------------------------------------------

void PDParticles::writeVTKFile(string tagname, int tagnum)
{

    // -----------------------------------
    //	Define the file location and name:
    // -----------------------------------

    ofstream outfile;
    std::stringstream filenamecombine;
    filenamecombine << "vtkoutput/" << tagname << "_" << tagnum << ".vtk";
    string filename = filenamecombine.str();
    outfile.open(filename.c_str(), ios::out | ios::app);

    // -----------------------------------
    //	Write the 'vtk' file header:
    // -----------------------------------

    string d = "   ";
    outfile << "# vtk DataFile Version 3.1" << endl;
    outfile << "VTK file containing particle data" << endl;
    outfile << "ASCII" << endl;
    outfile << " " << endl;
    outfile << "DATASET POLYDATA" << endl;
    outfile << " " << endl;
    outfile << "POINTS" << d << N << d << " float" << endl;

    // -----------------------------------
    //	Write the position data:
    // -----------------------------------

    for (int i=0; i<N; i++) {
        outfile << fixed << setprecision(3) << r[i*3+0] << d << r[i*3+1] << d << r[i*3+2] << endl;
    }

    // -----------------------------------
    //	write the radius data
    // -----------------------------------

    outfile << "POINT_DATA\t" << d << N << endl;
    outfile << "SCALARS radius float\n";
    outfile << "LOOKUP_TABLE default\n";

    for (int i=0; i<N; i++) {
        outfile << fixed << setprecision(3) << rad[i] << endl;
    }

    // -----------------------------------
    //	Close the file:
    // -----------------------------------

    outfile.close();
}



// -------------------------------------------------------------------------
// Setup the linked-list cell structure.
// -------------------------------------------------------------------------

void PDParticles::setupParticleCells()
{

   //	---------------------------------------
   // First, how many cells should exist:
   //	---------------------------------------

   cellWidth = rcut;
   ncellx = int(floor(box[0]/cellWidth));
   ncelly = int(floor(box[1]/cellWidth));
   ncellz = int(floor(box[2]/cellWidth));
   if (flag2D) ncellz = 1;

   // make sure domain is big enough
   if(!flag2D && (ncellz < 2 || ncelly < 2 || ncellx < 2))
   {
       cout << endl;
       if (ncellz < 2)
           cout << "\nDomain is not big enough in z-dim to support cell lists.\n";
       if (ncelly < 2)
           cout << "\nDomain is not big enough in y-dim to support cell lists.\n";
       if (ncellx < 2)
           cout << "\nDomain is not big enough in x-dim to support cell lists.\n";
       cout << "\nEither change the domain size or the cut off radius.\n";
       cout << endl << endl;
       throw 88;
   }
   else if (ncelly < 2 || ncellx < 2)
   {
       cout << endl;
       if (ncelly < 2)
           cout << "\nDomain is not big enough in y-dim to support cell lists.\n";
       if (ncellx < 2)
           cout << "\nDomain is not big enough in x-dim to support cell lists.\n";
       cout << "\nEither change the domain size or the cut off radius.\n";
       cout << endl << endl;
       throw 88;
   }

   ncell = ncellx*ncelly*ncellz;
   cellWidthx = box[0]/double(ncellx);
   cellWidthy = box[1]/double(ncelly);
   cellWidthz = box[2]/double(ncellz);

   //	---------------------------------------
   // Second, establish vector dimensions:
   //	---------------------------------------

   for (int i=0; i<ncell; i++) {
      head.push_back(-1);
      if (flag2D)  nncells = 4;
      if (!flag2D) nncells = 13;
      for (int i=0; i<nncells; i++) {
         cellmap.push_back(0);
      }
   }

   for (int i=0; i<N; i++) {
      list.push_back(0);
   }

   //	---------------------------------------
   // Last, create the cell nabor map:
   //	---------------------------------------

   for (int i=0; i<ncellx; i++) {
      for (int j=0; j<ncelly; j++) {
         for (int k=0; k<ncellz; k++) {

            int imap = cellIndex(i,j,k)*nncells;
            cellmap[imap+0]  = cellIndex(i+1, j , k );
            cellmap[imap+1]  = cellIndex(i+1,j+1, k );
            cellmap[imap+2]  = cellIndex( i ,j+1, k );
            cellmap[imap+3]  = cellIndex(i-1,j+1, k );
            if (nncells == 4) continue; // only do below if 3D cell structure
            cellmap[imap+4]  = cellIndex(i+1, j ,k-1);
            cellmap[imap+5]  = cellIndex(i+1,j+1,k-1);
            cellmap[imap+6]  = cellIndex( i ,j+1,k-1);
            cellmap[imap+7]  = cellIndex(i-1,j+1,k-1);
            cellmap[imap+8]  = cellIndex(i+1, j ,k+1);
            cellmap[imap+9]  = cellIndex(i+1,j+1,k+1);
            cellmap[imap+10] = cellIndex( i ,j+1,k+1);
            cellmap[imap+11] = cellIndex(i-1,j+1,k+1);
            cellmap[imap+12] = cellIndex( i , j ,k+1);

         }
      }
   }

}



// -------------------------------------------------------------------------
// Function to return the linked-list cell index given the
// x- y- z-coordinates of cell:
// -------------------------------------------------------------------------

int PDParticles::cellIndex(int i, int j, int k)
{
   if (i < 0) i += ncellx;
   if (i >= ncellx) i -= ncellx;
   if (j < 0) j += ncelly;
   if (j >= ncelly) j -= ncelly;
   if (k < 0) k += ncellz;
   if (k >= ncellz) k -= ncellz;
   return i*ncellz*ncelly + j*ncellz + k;
}



// -------------------------------------------------------------------------
// Function that calculates the total kinetic energy of the particles.
// -------------------------------------------------------------------------

double PDParticles::calcTotalKinEnergy()
{
    double ke = 0.0;
    double v2;
    for (int i=0; i<N;i++)
    {
        v2 = v[i*3]*v[i*3]+v[i*3+1]*v[i*3+1]+v[i*3+2]*v[i*3+2];
        ke += 0.5*mass[i]*v2;
    }
    return ke;
}



// -------------------------------------------------------------------------
// Function that writes the total kinetic energy of the particles to a
// file.
// -------------------------------------------------------------------------

void PDParticles::writeKinEnergy(std::vector<int> steps,std::vector<double> kinEnergy)
{
    ofstream ke;
    ke.open("equilibration.csv");
    // write header
    ke << "step,energy\n";
    for (size_t i = 0; i<steps.size();i++)
        ke << steps[i] << "," << kinEnergy[i] << "\n";
    ke.close();
}



// -------------------------------------------------------------------------
// Function that writes the particle forces and magnitudes to a csv file.
// -------------------------------------------------------------------------

void PDParticles::writeAllForces()
{
    string fname = "netForce";
    writeForce(current_step,f,fname);
}



// -------------------------------------------------------------------------
// Function that writes the particle forces and magnitudes to a csv file.
// -------------------------------------------------------------------------

void PDParticles::writeForce(int step,std::vector<double>& fr,string name)
{

    // -----------------------------------
    //	Define the file location and name:
    // -----------------------------------

    ofstream outfile;
    std::stringstream filenamecombine;
    filenamecombine << "vtkoutput/" << name << "_" << step << ".csv";
    string filename = filenamecombine.str();
    outfile.open(filename.c_str());

    // -----------------------------------
    //	write the file header
    // -----------------------------------

    outfile << "x,y,z,fx,fy,fz,mag\n";

    // -----------------------------------
    //	write force data
    // -----------------------------------

    double forceMag=0.0;
    for (int i = 0; i < N; i++)
    {
        forceMag = sqrt(fr[i*3]*fr[i*3]+fr[i*3+1]*fr[i*3+1]+fr[i*3+2]*fr[i*3+2]);
        outfile << r[i*3] << "," << r[i*3+1] << "," << r[i*3+2];
        outfile << "," << fr[i*3] << "," << fr[i*3+1] << "," << fr[i*3+2];
        outfile << "," << forceMag << "\n";
    }
}
