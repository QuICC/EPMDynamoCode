/** \file prec_FFTransform.cpp
 *  \brief Precision test for Fast Fourier transforms
 *
 *  \epmTSTodo Needs to be tested for MPI version
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
#include "GeneralScalars/RTPScalar.hpp"
#include "GeneralScalars/FFTFlatScalar.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Transforms/FFT/FFTransform.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;
/// Notation simplification typedef for FFT type
typedef epm::FFTransform  FFTransform;

/// Counter for the number of performed tests
int PERFORMED_TESTS = 0;

/// Error threshold value
epm::EPMFloat  ERROR_THRESHOLD = 1e-14;

/**
 * @brief Create test run values for forward/backward FFT transform loop test
 */
void initLoopTest(epm::FFTFlatScalar &rFFTValues, SmartTruncation pTrunc)
{
   // Make sure everything is 0
   rFFTValues.rFlat().setZero();

   // Fill used part with random values between (-1, 1)
   rFFTValues.rFlat().corner(Eigen::TopLeft, pTrunc->sim()->hoz()->nM(), rFFTValues.nCols()).setRandom();

   // Rescale values to (-1, 1)
   rFFTValues.rFlat().corner(Eigen::TopLeft, pTrunc->sim()->hoz()->nM(), rFFTValues.nCols()).array() -= 1.0;
   rFFTValues.rFlat().corner(Eigen::TopLeft, pTrunc->sim()->hoz()->nM(), rFFTValues.nCols()) *= 2.0;

   // Set m = 0 imaginary part to zero
   for(int c = 0; c < rFFTValues.nCols(); ++c)
   {
      rFFTValues.rFlat()(0, c).imag() = 0.0;
   }
}

/**
 * @brief Compare the two fft values
 */
int compareFFTValues(const epm::FFTFlatScalar &fftValues, const epm::FFTFlatScalar &fftValues2)
{
   // initialise failed comparisons counter
   int failed = 0;

   epm::EPMFloat  error;

   // Loop over all coefficients to count the wrong results
   for(int c = 0; c < fftValues.nCols(); ++c)
   {
      for(int m = 0; m < fftValues.nM() ; ++m)
      {
         // Compute difference between the real values
         error = std::abs(fftValues.flat()(m,c).real() - fftValues2.flat()(m,c).real());

         if(error > ERROR_THRESHOLD)
         {
            failed++;
         }

         // Compute difference between the imaginary values
         error = std::abs(fftValues.flat()(m,c).imag() - fftValues2.flat()(m,c).imag());

         if(error > ERROR_THRESHOLD)
         {
            failed++;
         }
      }
   }

   return failed;
}

/**
 * @brief Compute a tranform loop test
 */
int runLoopTest(SmartTruncation pTrunc, FFTransform &fft)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "FFT transform loop test" << std::endl;
   }

   // Initialise failed comparison counter
   int failed = 0;

   // Create RTPScalars
   epm::RTPScalar    rtpValues(pTrunc);

   // Create FFTFlatScalar
   epm::FFTFlatScalar   fftValues(pTrunc);
   epm::FFTFlatScalar   fftValues2(pTrunc);

   // Initialise the backward transform test values
   initLoopTest(fftValues, pTrunc); 

   // Compute backward transform
   fft.backward(rtpValues, fftValues);

   // Compute forward transform
   fft.forward(fftValues2, rtpValues);

   // Do the zero padding
   fftValues2.doZeroPadding();

   // Check result
   failed = compareFFTValues(fftValues, fftValues2);

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

   // Create the FFT object
   FFTransform    fft(pTrunc);

   int failed = 0;

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run forward FFT transform test
   failed += runLoopTest(pTrunc, fft);

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
