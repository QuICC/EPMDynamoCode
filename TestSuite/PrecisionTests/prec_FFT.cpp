/** \file prec_FFT.cpp
 *  \brief Precision test for Fast Fourier transforms
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

/**
 * @brief Create test run values for forward FFT transform test
 */
void initForwardTest(epm::RTPScalar &rRTPValues, epm::FFTFlatScalar &rFFTCorrect, SmartTruncation pTrunc)
{
   int m = 0;
   for(int r = 0; r < pTrunc->sim()->rad()->nR(); ++r)
   {
      for(int th = 0; th < pTrunc->sim()->hoz()->nTh(); ++th)
      {
         rRTPValues.rShell(r).col(th) = (static_cast<double>(m)*pTrunc->sim()->hoz()->phGrid()).cwise().cos() + (static_cast<double>(-m)*pTrunc->sim()->hoz()->phGrid()).cwise().sin() ;
         if(m < pTrunc->sim()->hoz()->nM())
         {
            m++;
         }
      }
   }

   int col = 0;
   m = 0;
   for(int r = 0; r < pTrunc->sim()->rad()->nR(); ++r)
   {
      for(int th = 0; th < pTrunc->sim()->hoz()->nTh(); ++th)
      {
         if(m != 0)
         {
            rFFTCorrect.rFlat()(m, col) = epm::EPMComplex(0.5, 0.5);
         } else
         {
            rFFTCorrect.rFlat()(m, col) = epm::EPMComplex(1.0, 0.0);
         }
         col++;
         if(m < pTrunc->sim()->hoz()->nM())
         {
            m++;
         }
      }
   }
}

/**
 * @brief Create test run values for forward FFT transform test
 */
void initBackwardTest(epm::FFTFlatScalar &rFFTValues, epm::RTPScalar &rRTPCorrect, SmartTruncation pTrunc)
{
   for(int m = 0; m < pTrunc->sim()->hoz()->nM(); ++m)
   {
      rFFTValues.rFlat().col(m).setRandom();
   }
}

/**
 * @brief Compute a forward FFT test transform
 */
int runForwardTest(SmartTruncation pTrunc, FFTransform &fft)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Initialise status counter
   int status = 0;

   // Create RTPScalar
   epm::RTPScalar    rtpValues(pTrunc);

   // Create FFTFlatScalars
   epm::FFTFlatScalar   fftValues(pTrunc);
   epm::FFTFlatScalar   fftCorrect(pTrunc);

   // Initialise the forward transform test values
   initForwardTest(rtpValues, fftCorrect, pTrunc); 

   // Compute forward transform
   fft.forward(fftValues, rtpValues);

   // Check result
   epm::EPMFloat  error;
   for(int c = 0; c < fftValues.nCols(); ++c)
   {
      error = (fftValues.flat().col(c) - fftCorrect.flat().col(c)).cwise().abs().maxCoeff();

      if(error > 1e-14)
      {
         status = 1;
      }
   }
   
   return status;
}

/**
 * @brief Compute a backward FFT test transform
 */
int runBackwardTest(SmartTruncation pTrunc, FFTransform &fft)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   // Initialise status counter
   int status = 1;

   // Create RTPScalars
   epm::RTPScalar    rtpValues(pTrunc);
   epm::RTPScalar    rtpCorrect(pTrunc);

   // Create FFTFlatScalar
   epm::FFTFlatScalar   fftValues(pTrunc);

   // Initialise the backward transform test values
   initBackwardTest(fftValues, rtpCorrect, pTrunc); 

   // Compute backward transform
   fft.backward(rtpValues, fftValues);

   // Check result

   return status;
}

/**
 * @brief Does transform test computation
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

   // Create the FFT object
   FFTransform    fft(pTrunc);

   int status = 0;

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run forward FFT transform test
   status += runForwardTest(pTrunc, fft);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run backward FFT transform test
   status += runBackwardTest(pTrunc, fft);

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
