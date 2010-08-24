/** \file unit_Truncation.cpp
 *  \brief Unit test for the truncation setup
 *
 *  \epmTSTodo Unit test not implemented yet
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
#include "Simulations/Types/WSHSimulation.hpp"
#include "Domain/Truncation.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;

/**
 * @brief Run the truncation unit test
 */
int runUnitTest()
{
   // Set test truncation values
   int maxN = 128;
   int maxL = 255;
   int maxM = 255;
   int Mp = 1;
   int nCore = 1444;

   // Initialise Spectral truncation information
   SmartTruncation   pTrunc = SimulationType::createTrunc(maxN, maxL, maxM, Mp, nCore);

   return 0;
}

/**
 * @brief General main, setting up MPI if required
 *
 * The actual program is in runUnitTest to make sure MPI initialisations
 * are called before anything else end finalization after destruction
 */
int main(int argc, char* argv[])
{
   // Initilise everything that can't be done inside a class
   epm::EPMDYNAMO_INITIALISER;

   // Storage for the return code of unit test
   int code;

   // Perform the actual test
   code = runUnitTest();

   // Finalise everything that can't be done inside a class
   epm::EPMDYNAMO_FINALIZER;

   // Create output to be taken up by CTest
   if(code)
   {
      std::cout << "Failed!" << std::endl;
   } else
   {
      std::cout << "Passed!" << std::endl;
   }

   return code;
}
