/**************************************************************************************/
/* 	                                                                                */
/*		MESO Code                                                                       */
/*		Authors: Prof. Paul Millett & members of the Millett Research Group             */
/* 	                                                                                */
/* 	Affiliation: University of Arkansas, Dept. of Mechanical Engineering            */
/*		Contact: pmillett@uark.edu, 479-575-2473                                        */
/* 	                                                                                */
/**************************************************************************************/

# include "MesoExecute.hpp"

int main(int argc, char **argv) {

//	----------------------------------------------------------
//	Create object that executes MESO simulation:
//	----------------------------------------------------------

	MesoExecute currentJob;

//	----------------------------------------------------------
//	Create simulation objects:
//	----------------------------------------------------------

	currentJob.createMesoObjects();

//	----------------------------------------------------------
//	Execute the simulation:
//	----------------------------------------------------------

	currentJob.executeMesoSimulation();

}
