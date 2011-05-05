/** \file SimulationRun.cpp
 *  \brief General simulation executable. Depending on preprocessor flag a different simulation will be built.
 */
#define SIMPATH Simulations/Implementations/IMPLEMENTATION.hpp
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )
#define SIMULATIONHEADER MAKE_STR( SIMPATH )
#define SIMULATION epm::IMPLEMENTATION

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iostream>

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Simulations/Simulation.hpp"
#include SIMULATIONHEADER
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"

namespace epm = EPMDynamo;

/**
 * @brief General setup for a simulation
 */
int runProgram()
{
   // Create the simulation object
   epm::Simulation<epm::WSHSimulation, SIMULATION>   sim;

   // Exception handling during the initialisation part
   try
   {
      // Initialise the whole simulation
      sim.init();
   }
   // If exception is thrown, finalise (close files) and return
   catch(int i)
   {
      // Cleanup and close file handles
      sim.finalise();

      return i;
   }

   // Run the simulation
   sim.run();

   // Cleanup and close file handles
   sim.finalise();

   return 0;
}

/**
 * @brief Main
 *
 * The actual program is in runProgram to make sure MPI initialisations
 * are called before anything else end finalization after destruction
 */
int main(int argc, char* argv[])
{
   // Initilise everything that can't be done inside a class
   epm::EPMDYNAMO_INITIALISER;

   int code;
   code = runProgram();

   // Finalise everything that can't be done inside a class
   epm::EPMDYNAMO_FINALIZER;

   return code;
}
