/** \file unit_SourceFile.cpp
 *  \brief Unit test for the reading/writing of HDF5 source files
 *
 *  \epmBug Unit test not implemented yet
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

namespace epm = EPMDynamo;

/**
 * @brief Run the HDF5 source file unit test
 */
int runUnitTest()
{
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
