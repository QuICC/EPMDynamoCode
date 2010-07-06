/** \file para_TubularSplitting.cpp
 *  \brief Parallelisation test for the tubular splitting algorithm
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
#include "General/EPMException.hpp"
#include "Domain/Truncation.hpp"
#include "Parallelisers/SplittingAlgorithms/TubularSplitting.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;

/**
 * @brief Run the truncation unit test
 */
int runParaTest()
{
   int maxN = 12;
   int maxL = 10;
   int maxM = 10;
   int Mp = 1;
   int nCore = 100;

   epm::SmartTruncation     pTrunc = SimulationType::createTrunc(maxN, maxL, maxM, Mp, 1);

   std::vector<epm::TubularSplitting>  pSplits;   


   for(int i = 0; i < nCore;  ++i)
   {
      pSplits.push_back(epm::TubularSplitting(pTrunc->sim(), nCore, i));
   }

   // RTP data storage
   int rR0;
   int rNr;
   epm::ArrayI rTh0;
   epm::ArrayI rNth;

   // FDSH data storage
   epm::ArrayI fR0;
   epm::ArrayI fNr;
   epm::ArrayI fM;

   // Spec data storage
   epm::ArrayI sL;
   std::vector<epm::ArrayI> sMs;

   for(int i = 0; i < nCore;  ++i)
   {
      // Compute RTP data box splitting
      pSplits.at(i).splitRTP(rR0, rNr, rTh0, rNth, i);

      std::cerr << "r0: " << rR0 << std::endl;
      std::cerr << "nR: " << rNr << std::endl;
      std::cerr << "th0: " << rTh0.transpose() << std::endl;
      std::cerr << "nTh: " << rNth.transpose() << std::endl;

      // Compute FDSH data box splitting
      pSplits.at(i).splitFDSH(fR0, fNr, fM, i);
      std::cerr << "fR0: " << fR0.transpose() << std::endl;
      std::cerr << "fNr: " << fNr.transpose() << std::endl;
      std::cerr << "fM: " << fM.transpose() << std::endl;

      // Compute Spectral data box splitting
      pSplits.at(i).splitSpec(sL, sMs, i);
      std::cerr << "sL: " << sL.transpose() << std::endl;
      for(unsigned int j = 0; j < sMs.size(); ++j)
      {
         std::cerr << "sMs(" << j << "): " << sMs.at(j).transpose() << std::endl;
      }
   }

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
   int code = 1;

   try
   {
      // Perform the actual test
      code = runParaTest();
   } 
   catch(epm::EPMException &e)
   {
      e.printStdMessage();
   }

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
