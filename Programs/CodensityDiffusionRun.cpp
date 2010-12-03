/** \file CodensityDiffusionRun.cpp
 *  \brief Codensity diffusion simulation executable
 */

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iostream>

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Simulations/Simulation.hpp"
#include "Simulations/Implementations/CodensityDiffusionSimulation.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/Types/WSHSimInc.hpp"

namespace epm = EPMDynamo;

/**
 * @brief Codensity diffusion simulation
 */
int runProgram()
{
   // Create the simulation object
   epm::Simulation<epm::WSHSimulation, epm::CodensityDiffusionSimulation>   sim;

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
