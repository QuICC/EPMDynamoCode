/** \file unit_Worland.cpp
 *  \brief Unit test for the Worland polynomials
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
#include "General/EPMException.hpp"
#include "Domain/Truncation.hpp"
#include "Polynomials/RadialBasis.hpp"
#include "Polynomials/WorlandPolynomial.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;
/// Notation simplification typedef for radial basis type
typedef SimulationType::RadialBasisType RadialBasisType;

/**
 * @brief Run the Worland unit test
 */
int runUnitTest()
{
   // Set some truncation values
   int maxN = 12;
   int maxL = 32;
   int maxM = 10;
   int Mp = 1;
   int nCore = 1;

   const double ERROR_THRESHOLD = 1e-13;

   SmartTruncation  pTrunc = SimulationType::createTrunc(maxN, maxL, maxM, Mp, nCore);

   int nN = pTrunc->sim()->rad()->nN();
   int nL = pTrunc->local()->spec()->nL();

   RadialBasisType radBasis(pTrunc);

   epm::EPMFloat tmp = 0.0;
   epm::EPMFloat maxCoeff = 0.0;

   for(int l = 0; l < nL ; ++l)
   {
      tmp = ((radBasis.at(l).poly() * radBasis.at(l).wPoly()) - epm::Matrix::Identity(nN, nN)).maxCoeff();
      maxCoeff = std::max(tmp, maxCoeff);
   }

   if(maxCoeff < ERROR_THRESHOLD)
   {
      return 0;
   } else
   {
      return 1;
   }
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
   int code = -1;

   try
   {
      // Perform the actual test
      code = runUnitTest();
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
