/** \file prec_LegendreTransform.cpp
 *  \brief Precision test for Legendre transforms
 *
 *  \epmTSTodo Needs to be implemented correctly
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
#include "GeneralScalars/FDSHOrderScalar.hpp"
#include "GeneralScalars/FFTOrderScalar.hpp"
#include "Polynomials/AssociatedLegendreOperator.hpp"
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

/// Error threshold value
epm::EPMFloat  ERROR_THRESHOLD = 1e-10;

//
// Loop test initialisation
//

/**
 * @brief Initialise test values for transform loop test
 */
void initLoopTest(epm::FDSHOrderScalar &rFDSHValues, SmartTruncation pTrunc)
{
   // Loop over harmonic orders shells
   for(int m = 0; m < rFDSHValues.nM(); ++m)
   {
      // Set random values
      rFDSHValues.rMShell(m).setRandom();

      // Rescale values to (-1 ,1) interval
      rFDSHValues.rMShell(m).cwise() -= 1.0;
      rFDSHValues.rMShell(m) *= 2.0;
   }

   // Set imaginary part of m=0 to zero
   for(int r = 0; r < rFDSHValues.nR(0); ++r)
   {
      for(int l = 0; l < rFDSHValues.nL(0); ++l)
      {
         rFDSHValues.rMShell(0)(l, r).imag() = 0.0;
      }
   }
}


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
// Test transform loop
//

/**
 * @brief Compare the two fdsh values
 */
int compareFDSHValues(const epm::FDSHOrderScalar &fdshValues, const epm::FDSHOrderScalar &fdshValues2)
{
   // initialise failed comparisons counter
   int failed = 0;

   epm::EPMFloat  error;

   // Loop over all coefficients to count the wrong results
   for(int m = 0; m < fdshValues.nM() ; ++m)
   {
      for(int r = 0; r < fdshValues.nR(m); ++r)
      {
         for(int l = 0; l < fdshValues.nL(m); ++l)
         {
            // Compute difference between real values
            error = std::abs(fdshValues.mshell(m)(l,r).real() - fdshValues2.mshell(m)(l,r).real());

            if(error > ERROR_THRESHOLD)
            {
               failed++;
            }

            // Compute difference between real values
            error = std::abs(fdshValues.mshell(m)(l,r).imag() - fdshValues2.mshell(m)(l,r).imag());

            if(error > ERROR_THRESHOLD)
            {
               failed++;
            }
         }
      }
   }

   return failed;
}

/**
 * @brief Test run for the transform loop
 */
int runLoopTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Legendre transform loop test" << std::endl;
   }

   int failed = 0;

   // Create FFTOrderScalar
   epm::FFTOrderScalar  fftValues(pTrunc);

   // Create FDSHOrderScalars
   epm::FDSHOrderScalar fdshValues(pTrunc);
   epm::FDSHOrderScalar fdshValues2(pTrunc);

   // Initialise the forward transform test values
   initLoopTest(fdshValues, pTrunc); 

   // Compute the projection
   legTransform.multM<epm::SetProduct>(fftValues, fdshValues, &epm::AssociatedLegendreOperator::proj);
   
   // Compute the integration
   legTransform.multM<epm::SetProduct>(fdshValues2, fftValues, &epm::AssociatedLegendreOperator::intg);

   // Check result
   failed = compareFDSHValues(fdshValues, fdshValues2);

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
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

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intg\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgTest(pTrunc); 

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "intgTh2S" integration
 */
int runIntgTh2STest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgTh2S\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgTh2STest(pTrunc); 

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "intgTh2T" integration
 */
int runIntgTh2TTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgTh2T\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgTh2TTest(pTrunc); 

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "intgPh2S" integration
 */
int runIntgPh2STest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgPh2S\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgPh2STest(pTrunc); 

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "intgPh2T" integration
 */
int runIntgPh2TTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgPh2T\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgPh2TTest(pTrunc); 

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
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

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"proj\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjTest(pTrunc); 

   // Compute transform

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "projS2Th" projection
 */
int runProjS2ThTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projS2Th\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjS2ThTest(pTrunc); 

   // Compute transform

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "projT2Th" projection
 */
int runProjT2ThTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projT2Th\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjT2ThTest(pTrunc); 

   // Compute transform

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "projS2Ph" projection
 */
int runProjS2PhTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projS2Ph\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjS2PhTest(pTrunc); 

   // Compute transform

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "projT2Ph" projection
 */
int runProjT2PhTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projT2Ph\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjT2PhTest(pTrunc); 

   // Compute transform

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "projGradPh" projection
 */
int runProj2GradPhTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"proj2GradPh\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProj2GradPhTest(pTrunc); 

   // Compute transform

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Test run for the "proj2GradPh" projection
 */
int runProj2GradThTest(SmartTruncation pTrunc, LegendreTransform &legTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"proj2GradTh\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProj2GradPhTest(pTrunc); 

   // Compute transform

   // Check result

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}


/**
 * @brief Perform all the test computations
 */
int runPrecTest()
{
   // Set test truncation values
   int maxN = 100;
   int maxL = 256;
   int maxM = 256;
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

   int failed = 0;

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intg" test
   failed += runLoopTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intg" test
   failed += runIntgTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgTh2S" test
   failed += runIntgTh2STest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgTh2T" test
   failed += runIntgTh2TTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgPh2S" test
   failed += runIntgPh2STest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgPh2T" test
   failed += runIntgPh2TTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj" test
   failed += runProjTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projS2Th" test
   failed += runProjS2ThTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projT2Th" test
   failed += runProjT2ThTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projS2Ph" test
   failed += runProjS2PhTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projT2Ph" test
   failed += runProjT2PhTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj2GradPh" test
   failed += runProj2GradPhTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj2GradTh" test
   failed += runProj2GradThTest(pTrunc, legTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   return failed;
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
