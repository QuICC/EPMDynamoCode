/** \file prec_WorlandTransform.cpp
 *  \brief Precision test for Worland transforms
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

// TestSuite includes
//
#include "ErrorComputer.hpp"

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "General/MathConstants.hpp"
#include "Domain/Truncation.hpp"
#include "Polynomials/WorlandPolynomial.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Transforms/Polynomial/RadialTransform.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;
/// Notation simplification typedef for Worland transform type
typedef epm::RadialTransform<epm::WorlandPolynomial> WorlandTransform;

/// Counter for the number of performed tests
int PERFORMED_TESTS = 0;

/// Error threshold value
epm::EPMFloat  ERROR_THRESHOLD = 1e-12;

//
// Loop test initialisation
//

/**
 * @brief Initialise test values for transform loop test
 */
void initLoopTest(epm::SpectralSHScalar &rSSHValues, SmartTruncation pTrunc)
{
   // Loop over harmonic degree shells
   for(int l = 0; l < rSSHValues.nL(); ++l)
   {
      // Set random values
      rSSHValues.rLShell(l).setRandom();

      // Rescale values to (-1 ,1) interval
      rSSHValues.rLShell(l).array() -= 1.0;
      rSSHValues.rLShell(l) *= 2.0;
   }

   // Set imaginary part of m=0 to zero
   for(int l = 0; l < rSSHValues.nL(); ++l)
   {
      for(int n = 0; n < rSSHValues.nN(); ++n)
      {
         rSSHValues.rLShell(l)(n, 0).imag() = 0.0;
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
 * @brief Initialise test values for "intgQ2Pol"
 */
void initIntgQ2PolTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgT2Tor"
 */
void initIntgT2TorTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgQ2CurlCurlProj"
 */
void initIntgQ2CurlCurlProjTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgS2CurlCurlProj"
 */
void initIntgS2CurlCurlProjTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "intgT2CurlProj"
 */
void initIntgT2CurlProjTest(SmartTruncation pTrunc)
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
 * @brief Initialise test values for "proj2GradTP"
 */
void initProj2GradTPTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projPol2S"
 */
void initProjPol2STest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projTor2T"
 */
void initProjTor2TTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projTor2CurlQ"
 */
void initProjTor2CurlQTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projTor2CurlS"
 */
void initProjTor2CurlSTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projPol2CurlT"
 */
void initProjPol2CurlTTest(SmartTruncation pTrunc)
{
}

/**
 * @brief Initialise test values for "projPol2Q"
 */
void initProjPol2QTest(SmartTruncation pTrunc)
{
}

//
// Test transform loop
//

/**
 * @brief Compare the two fdsh values
 */
int compareSSHValues(const epm::SpectralSHScalar &sshValues, const epm::SpectralSHScalar &sshValues2)
{
   // initialise failed comparisons counter
   int failed = 0;

   epm::EPMFloat  error;

   // Loop over all coefficients to count the wrong results
   for(int l = 0; l < sshValues.nL() ; ++l)
   {
      for(int m = 0; m < sshValues.nM(l); ++m)
      {
         for(int n = 0; n < sshValues.nN(); ++n)
         {
            // Compute difference between real values
            error = std::abs(sshValues.lshell(l)(n,m).real() - sshValues2.lshell(l)(n,m).real());

            if(error > ERROR_THRESHOLD)
            {
               failed++;
            }

            // Compute difference between real values
            error = std::abs(sshValues.lshell(l)(n,m).imag() - sshValues2.lshell(l)(n,m).imag());

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
int runLoopTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Worland transform loop test" << std::endl;
   }

   int failed = 0;

   // Create FDSHDegreeScalar
   epm::FDSHDegreeScalar  fdshValues(pTrunc);

   // Create SpectralSHScalars
   epm::SpectralSHScalar sshValues(pTrunc, true);
   epm::SpectralSHScalar sshValues2(pTrunc, true);

   // Initialise the forward transform test values
   initLoopTest(sshValues, pTrunc); 

   // Compute the projection
   wTransform.multL<epm::SetProduct>(fdshValues, sshValues, &epm::TorPolRadialOperator<epm::WorlandPolynomial>::proj);
   
   // Compute the integration
   wTransform.multL<epm::SetProduct>(sshValues2, fdshValues, &epm::TorPolRadialOperator<epm::WorlandPolynomial>::intg);

   // Check result
   failed = compareSSHValues(sshValues, sshValues2);

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
int runIntgTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
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
 * @brief Test run for the "intgQ2Pol" integration
 */
int runIntgQ2PolTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgQ2Pol\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgQ2PolTest(pTrunc); 

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
 * @brief Test run for the "intgT2Tor" integration
 */
int runIntgT2TorTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgT2Tor\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgT2TorTest(pTrunc); 

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
 * @brief Test run for the "intgQ2CurlCurlProj" integration
 */
int runIntgQ2CurlCurlProjTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgQ2CurlCurlProj\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgQ2CurlCurlProjTest(pTrunc); 

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
 * @brief Test run for the "intgS2CurlCurlProj" integration
 */
int runIntgS2CurlCurlProjTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgS2CurlCurlProj\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgS2CurlCurlProjTest(pTrunc); 

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
 * @brief Test run for the "intgT2CurlProj" integration
 */
int runIntgT2CurlProjTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"intgT2CurlProj\" integration test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initIntgT2CurlProjTest(pTrunc); 

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


//
// Test runs for projection computations
//

/**
 * @brief Test run for the "proj" projection
 */
int runProjTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
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
 * @brief Test run for the "proj2GradTP" projection
 */
int runProj2GradTPTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"proj2GradTP\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProj2GradTPTest(pTrunc); 

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
 * @brief Test run for the "projPol2S" projection
 */
int runProjPol2STest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projPol2S\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjPol2STest(pTrunc); 

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
 * @brief Test run for the "projTor2T" projection
 */
int runProjTor2TTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projTor2T\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjTor2TTest(pTrunc); 

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
 * @brief Test run for the "projTor2CurlQ" projection
 */
int runProjTor2CurlQTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projTor2CurlQ\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjTor2CurlQTest(pTrunc); 

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
 * @brief Test run for the "projTor2CurlS" projection
 */
int runProjTor2CurlSTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projTor2CurlS\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjTor2CurlSTest(pTrunc); 

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
 * @brief Test run for the "projPol2CurlT" projection
 */
int runProjPol2CurlTTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projPol2CurlT\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjPol2CurlTTest(pTrunc); 

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
 * @brief Test run for the "projPol2Q" projection
 */
int runProjPol2QTest(SmartTruncation pTrunc, WorlandTransform &wTransform)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "\"projPol2Q\" projection test" << std::endl;
   }

   int failed = 1;

   // Initialise the forward transform test values
   initProjPol2QTest(pTrunc); 

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
 * @brief Does transform test computation
 */
int runPrecTest()
{
   // Set test truncation values
   int maxN = 100;
   int maxL = 64;
   int maxM = 64;
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

   // Create the Worland transform object
   WorlandTransform    wTransform(pTrunc);

   int failed = 0;

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intg" test
   failed += runLoopTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intg" test
   failed += runIntgTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgQ2Pol" test
   failed += runIntgQ2PolTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgT2Tor" test
   failed += runIntgT2TorTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgQ2CurlCurlProj" test
   failed += runIntgQ2CurlCurlProjTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgS2CurlCurlProj" test
   failed += runIntgS2CurlCurlProjTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "intgT2CurlProj" test
   failed += runIntgT2CurlProjTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj" test
   failed += runProjTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "proj2GradTP" test
   failed += runProj2GradTPTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projPol2S" test
   failed += runProjPol2STest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projTor2T" test
   failed += runProjTor2TTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projTor2CurlQ" test
   failed += runProjTor2CurlQTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projTor2CurlS" test
   failed += runProjTor2CurlSTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projPol2CurlT" test
   failed += runProjPol2CurlTTest(pTrunc, wTransform);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run "projPol2Q" test
   failed += runProjPol2QTest(pTrunc, wTransform);

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
