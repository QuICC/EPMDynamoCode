/** \file prec_AssocLegendreBasis.cpp
 *  \brief Precision test for the Associated Legendre polynomials
 *
 *  \epmBug Needs to be implemented correctly
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
#include "Polynomials/AssocLegendreBasis.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;

/// Base name for all files
std::string testBasename = "alp_";

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

   // write Harmonic truncation to file
   file << pTrunc->sim()->hoz()->maxL() << std::endl;
   file << pTrunc->sim()->hoz()->maxM();

   // close file
   file.close();

   return status;
}

/**
 * @brief Write "grid" data to file
 */
int writeGridData(SmartTruncation pTrunc, epm::AssocLegendreBasis &legBasis)
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
   converter << pTrunc->sim()->hoz()->nTh();

   // Create base name
   basename = testBasename + "N" + converter.str() + "_grid.dat";

   // Empty string stream
   converter.str("");

   // Open file
   file.open(basename.c_str());

   // Write data to file
   file << legBasis.at(0).grid();

   // close file
   file.close();

   return status;
}

/**
 * @brief Write "weights" data to file
 */
int writeWeightsData(SmartTruncation pTrunc, epm::AssocLegendreBasis &legBasis)
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
   converter << pTrunc->sim()->hoz()->nTh();

   // Create base name
   basename = testBasename + "N" + converter.str() + "_weights.dat";

   // Empty string stream
   converter.str("");

   // Open file
   file.open(basename.c_str());

   // Write data to file
   file << legBasis.at(0).weights();

   // close file
   file.close();

   return status;
}

/**
 * @brief Write "intg" data to file
 */
int writeIntgData(SmartTruncation pTrunc, epm::AssocLegendreBasis &legBasis)
{
   // File writing status
   int status = 0;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Writing \"intg\" data to file" << std::endl;
   }

   // Create file object
   std::ofstream  file;

   // Filename base string
   std::string basename;

   // Filename string
   std::string filename;

   // Stringstream for the conversion
   std::stringstream converter;

   // Convert grid size
   converter << pTrunc->sim()->hoz()->nTh();

   // Create base name
   basename = testBasename + "N" + converter.str() + "_intg" + "_M";

   // Empty string stream
   converter.str("");

   // Loop over all harmonic orders
   for(int m = 0; m < pTrunc->sim()->hoz()->nM(); ++m)
   {
      // Convert harmonic order
      converter << m;

      // Build filename
      filename = basename + converter.str() + ".dat";

      // Open file
      file.open(filename.c_str());

      // Write data to file
      file << legBasis.at(m).intg().op();

      // close file
      file.close();

      // empty stringstream
      converter.str("");
   }

   return status;
}

/**
 * @brief Write "proj" data to file
 */
int writeProjData(SmartTruncation pTrunc, epm::AssocLegendreBasis &legBasis)
{
   // File writing status
   int status = 0;

   // Test presentation output
   if(pTrunc->para().id() == 0)
   {
      std::cout << "Writing \"proj\" data to file" << std::endl;
   }

   // Create file object
   std::ofstream  file;

   // Filename base string
   std::string basename;

   // Filename string
   std::string filename;

   // Stringstream for the conversion
   std::stringstream converter;

   // Convert grid size
   converter << pTrunc->sim()->hoz()->nTh();

   // Create base name
   basename = testBasename + "N" + converter.str() + "_proj" + "_M";

   // Empty string stream
   converter.str("");

   // Loop over all harmonic orders
   for(int m = 0; m < pTrunc->sim()->hoz()->nM(); ++m)
   {
      // Convert harmonic order
      converter << m;

      // Build filename
      filename = basename + converter.str() + ".dat";

      // Open file
      file.open(filename.c_str());

      // Write data to file
      file << legBasis.at(m).proj().op();

      // close file
      file.close();

      // empty stringstream
      converter.str("");
   }

   return status;
}

/**
 * @brief Perform all the tests
 */
int runPrecTest()
{
   // Set test truncation values
   int maxN = 100;
   int maxL = 64;
   int maxM = 64;
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
      std::cout << "\t Theta N: " << pTrunc->sim()->hoz()->nTh() << std::endl;
      std::cout << "\t Phi N: " << pTrunc->sim()->hoz()->nPh() << std::endl;
      std::cout << "\t Max N: " << pTrunc->sim()->rad()->maxN() << std::endl;
      std::cout << "\t Max L: " << pTrunc->sim()->hoz()->maxL() << std::endl;
      std::cout << "\t Max M: " << pTrunc->sim()->hoz()->maxM() << std::endl;
      std::cout << "\t Mp : " << pTrunc->sim()->hoz()->mp() << std::endl;
      std::cout << std::endl << std::endl;
   }

   // Create the Associated Legendre Basis object
   epm::AssocLegendreBasis    legBasis(pTrunc);

   // Write setup data to file
   status += writeSetupData(pTrunc);

   // Write "grid" data to file
   status += writeGridData(pTrunc, legBasis);

   // Write "weights" data to file
   status += writeWeightsData(pTrunc, legBasis);

   // Write "intg" data to file
   status += writeIntgData(pTrunc, legBasis);

   // Write "proj" data to file
   status += writeProjData(pTrunc, legBasis);

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
