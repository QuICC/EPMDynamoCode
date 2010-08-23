/** \file prec_WorlandBasis.cpp
 *  \brief Precision test for the Worland polynomials
 *
 *  \epmTSTodo Needs to be implemented correctly
 */

#define EIGEN_DEFAULT_IO_FORMAT IOFormat(16)

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iostream>
#include <fstream>
#include <sstream>
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
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/WorlandPolynomial.hpp"
#include "Polynomials/RadialBasis.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;
/// Notation simplification typedef for worland operator type
typedef epm::TorPolRadialOperator<epm::WorlandPolynomial> WorlandOperator;
/// Notation simplification typedef for worland radial basis type
typedef epm::RadialBasis<WorlandOperator> WorlandBasis;

/// Base name for all files
std::string testBasename = "wp_";

/// Base name for the grid data files
std::string gridBasename;

/// File ordering character
const std::string orderString = "_L";

/**
 * @brief Write setup data to file
 */
int writeSetupData(SmartTruncation pTrunc)
{
   // File writing status
   int status = 0;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Writing setup data to file" << std::endl;
   }

   // Create file object
   std::ofstream  file;

   // Filename base string
   std::string filename;

   // Create base name
   filename = testBasename + "setup.dat";

   // Open file
   file.open(filename.c_str());

   // write Harmonic and radial truncation to file
   file << pTrunc->sim()->rad()->maxN() << std::endl;
   file << pTrunc->sim()->hoz()->maxL();

   // close file
   file.close();

   return status;
}

/**
 * @brief Create grid data filename base
 */
void setGridBasename(SmartTruncation pTrunc)
{
   // Stringstream for the conversion
   std::stringstream converter;

   // Convert grid size
   converter << pTrunc->sim()->rad()->nR();

   // Create base name
   gridBasename = testBasename + "N" + converter.str() + "_";

   // Empty string stream
   converter.str("");
}

/**
 * @brief Write "grid" data to file
 */
int writeGridData(SmartTruncation pTrunc, WorlandBasis &wBasis)
{
   // File writing status
   int status = 0;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Writing \"grid\" data to file" << std::endl;
   }

   // Create file object
   std::ofstream  file;

   // Filename base string
   std::string basename;

   // Stringstream for the conversion
   std::stringstream converter;

   // Convert grid size
   converter << pTrunc->sim()->rad()->nR();

   // Create base name
   basename = gridBasename + "grid.dat";

   // Empty string stream
   converter.str("");

   // Open file
   file.open(basename.c_str());

   // Write data to file
   file << wBasis.at(0).grid();

   // close file
   file.close();

   return status;
}

/**
 * @brief Write data to file
 */
template<class T, typename TExp> int writeData(std::string name, const epm::PolynomialOperator<TExp>& (T::*Tptr)() const, SmartTruncation pTrunc, WorlandBasis &wBasis)
{
   // File writing status
   int status = 0;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Writing \"" + name + "\" data to file" << std::endl;
   }

   // Create file object
   std::ofstream  file;

   // Filename base string
   std::string basename;

   // Filename string
   std::string filename;

   // Stringstream for the conversion
   std::stringstream converter;

   // Create base name
   basename = gridBasename + name + orderString;

   // Empty string stream
   converter.str("");

   // Loop over all harmonic orders
   for(int l = 0; l < pTrunc->sim()->hoz()->nL(); ++l)
   {
      // Convert harmonic order
      converter << l;

      // Build filename
      filename = basename + converter.str() + ".dat";

      // Open file
      file.open(filename.c_str());

      // Write data to file
      file << (wBasis.at(l).*Tptr)().op();

      // close file
      file.close();

      // empty stringstream
      converter.str("");
   }

   return status;
}

/**
 * @brief Write "weights" data to file
 */
int writeWeightsData(SmartTruncation pTrunc, WorlandBasis &wBasis)
{
   // File writing status
   int status = 0;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Writing \"weights\" data to file" << std::endl;
   }

   // Create file object
   std::ofstream  file;

   // Filename base string
   std::string basename;

   // Stringstream for the conversion
   std::stringstream converter;

   // Convert grid size
   converter << pTrunc->sim()->rad()->nR();

   // Create base name
   basename = gridBasename + "weights.dat";

   // Empty string stream
   converter.str("");

   // Open file
   file.open(basename.c_str());

   // Write data to file
   file << wBasis.at(0).weights();

   // close file
   file.close();

   return status;
}

/**
 * @brief Perform all the tests
 */
int runPrecTest()
{
   // Set test truncation values
   int maxN = 20;
   int maxL = 32;
   int maxM = 32;
   int Mp = 1;
   int nCore = 1;

   int status = 0;

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
      std::cout << "\t Theta N: " << pTrunc->sim()->rad()->nR() << std::endl;
      std::cout << "\t Phi N: " << pTrunc->sim()->hoz()->nPh() << std::endl;
      std::cout << "\t Max N: " << pTrunc->sim()->rad()->maxN() << std::endl;
      std::cout << "\t Max L: " << pTrunc->sim()->hoz()->maxL() << std::endl;
      std::cout << "\t Max M: " << pTrunc->sim()->hoz()->maxM() << std::endl;
      std::cout << "\t Mp : " << pTrunc->sim()->hoz()->mp() << std::endl;
      std::cout << std::endl << std::endl;
   }

   // Create the Worland basis object
   WorlandBasis    wBasis(pTrunc);

   // Set "grid" data basename
   setGridBasename(pTrunc);

   // Write setup data to file
   status += writeSetupData(pTrunc);

   // Write "intg" data to file
   status += writeData("intg", &WorlandOperator::intg, pTrunc, wBasis);

   // Write "intgQ2Pol" data to file
   status += writeData("intgQ2Pol", &WorlandOperator::intgQ2Pol, pTrunc, wBasis);

   // Write "intgT2Tor" data to file
   status += writeData("intgT2Tor", &WorlandOperator::intgT2Tor, pTrunc, wBasis);

   // Write "intgQ2CurlCurlProj" data to file
   status += writeData("intgQ2CurlCurlProj", &WorlandOperator::intgQ2CurlCurlProj, pTrunc, wBasis);

   // Write "intgS2CurlCurlProj" data to file
   status += writeData("intgS2CurlCurlProj", &WorlandOperator::intgS2CurlCurlProj, pTrunc, wBasis);

   // Write "intgT2CurlProj" data to file
   status += writeData("intgT2CurlProj", &WorlandOperator::intgT2CurlProj, pTrunc, wBasis);

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
         std::cout << "Failed!" << std::endl;
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
