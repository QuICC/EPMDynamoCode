/** \file prec_LegendreTransform.cpp
 *  \brief Precision test for Legendre transforms
 *
 *  \epmBug Needs to be implemented correctly
 */

#define EIGEN_DEFAULT_IO_FORMAT IOFormat(10)

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iostream>
#include <Eigen/Array>

// TestSuite includes
//
#include "ErrorComputer.hpp"

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "General/MathConstants.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FFTFlatScalar.hpp"
#include "GeneralScalars/FFTOrderScalar.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Transforms/Polynomial/LegendreTransform.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;
/// Notation simplification typedef for Legendre transform type
typedef epm::LegendreTransform  LegendreTransform;

/// Counter for the number of performed tests
int PERFORMED_TESTS = 0;

//
// Integration tests initialisations
//

/**
 * @brief Initialise test values for "intg"
 */
void initIntgTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgTh2S"
 */
void initIntgTh2STest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgTh2T"
 */
void initIntgTh2TTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgPh2S"
 */
void initIntgPh2STest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgPh2T"
 */
void initIntgPh2TTest(SmartTruncation pTrunc)
{
}


//
// Projection tests initialisations
//

/**
 * @brief Initialise test values for "proj"
 */
void initProjTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projS2Th"
 */
void initProjS2ThTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projT2Th"
 */
void initProjT2ThTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projS2Ph"
 */
void initProjS2PhTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projT2Ph
 */
void initProjT2PhTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "proj2GradPh"
 */
void initProj2GradPhTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "proj2GradTh"
 */
void initProj2GradThTest(SmartTruncation pTrunc)
{
}


//
// Test runs for integration computations
//

/**
 * @brief Test run for the "intg" integration
 */
int runIntgTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initIntgTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "intgTh2S" integration
 */
int runIntgTh2STest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initIntgTh2STest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "intgTh2T" integration
 */
int runIntgTh2TTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initIntgTh2TTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "intgPh2S" integration
 */
int runIntgPh2STest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initIntgPh2STest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "intgPh2T" integration
 */
int runIntgPh2TTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initIntgPh2TTest(pTrunc); 

   return status;
}


//
// Test runs for projection computations
//

/**
 * @brief Test run for the "proj" projection
 */
int runProjTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initProjTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "projS2Th" projection
 */
int runProjS2ThTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initProjS2ThTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "projT2Th" projection
 */
int runProjT2ThTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initProjT2ThTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "projS2Ph" projection
 */
int runProjS2PhTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initProjS2PhTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "projT2Ph" projection
 */
int runProjT2PhTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initProjT2PhTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "projGradPh" projection
 */
int runProj2GradPhTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initProj2GradPhTest(pTrunc); 

   return status;
}

/**
 * @brief Test run for the "proj2GradPh" projection
 */
int runProj2GradThTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   int status = 1;

   // Initialise the forward transform test values
   initProj2GradPhTest(pTrunc); 

   return status;
}


/**
 * @brief Perform all the test computations
 */
int runPrecTest()
{
   // Set test truncation values
   int maxN = 20;
   int maxL = 32;
   int maxM = 32;
   int Mp = 1;
   int nCore = 1;

   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Comm_size(MPI_COMM_WORLD, &nCore);
   #endif //EPMDYNAMO_MPI

   // Initialise Spectral truncation information
   SmartTruncation   pTrunc = SimulationType::createTrunc(maxN, maxL, maxM, Mp, nCore);

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Truncation information:" << std::endl;
      std::cout << "\t Radial N: " << pTrunc->sim()->rad()->nR() << std::endl;
      std::cout << "\t Theta N: " << pTrunc->sim()->hoz()->nTh() << std::endl;
      std::cout << "\t Phi N: " << pTrunc->sim()->hoz()->nPh() << std::endl;
      std::cout << "\t Max N: " << pTrunc->sim()->rad()->maxN() << std::endl;
      std::cout << "\t Max L: " << pTrunc->sim()->hoz()->maxL() << std::endl;
      std::cout << "\t Max M: " << pTrunc->sim()->hoz()->maxM() << std::endl;
      std::cout << "\t Mp : " << pTrunc->sim()->hoz()->mp() << std::endl;
      std::cout << std::endl << std::endl;
   }

   // Create the Legendre transform object
   LegendreTransform    legTransform(pTrunc);

   int status = 0;

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intg" test
   status += runIntgTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgTh2S" test
   status += runIntgTh2STest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgTh2T" test
   status += runIntgTh2TTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgPh2S" test
   status += runIntgPh2STest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgPh2T" test
   status += runIntgPh2TTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj" test
   status += runProjTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projS2Th" test
   status += runProjS2ThTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projT2Th" test
   status += runProjT2ThTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projS2Ph" test
   status += runProjS2PhTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projT2Ph" test
   status += runProjT2PhTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj2GradPh" test
   status += runProj2GradPhTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj2GradTh" test
   status += runProj2GradThTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   return status;
}

/**
 * @brief General main, setting up MPI if required
 *
 * The actual program is in runPrecTest to make sure MPI initialisations
 * are called before anything else end finalization after destruction
 */
int main(int argc, char* argv[])
{
   // Initilise everything that can't be done inside a class
   epm::EPMDYNAMO_INITIALISER;

   int code = -1;

   try
   {
      code = runPrecTest();
   }
   catch(epm::EPMException &e)
   {
      e.printStdMessage();
   }

   // Create output to be taken up by CTest
   if(code)
   {
      if(epm::EPMDYNAMO_RANK == 0)
      {
         std::cout << "Failed! (" << code << " out of " << PERFORMED_TESTS << ")" << std::endl;
      }
   } else
   {
      if(epm::EPMDYNAMO_RANK == 0)
      {
         std::cout << "Passed!" << std::endl;
      }
   }

   // Finalise everything that can't be done inside a class
   epm::EPMDYNAMO_FINALIZER;

   return code;
}
