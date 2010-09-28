/** \file TransformOperatorTest.hpp
 *  \brief Simple base class for the transform operators used for the testsuite.
 */

#ifndef TRANSFORMOPERATORTEST_HPP
#define TRANSFORMOPERATORTEST_HPP

// System includes
//
#include <fstream>
#include <sstream>
#include <string>

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Polynomials/PolynomialOperator.hpp"

namespace EPMDynamo {

   /**
    * @brief Simple base class for the transform operators used for the testsuite.
    */
   template <typename TBasis> class TransformOperatorTest
   {
      public:
         /**
         * @brief Constructs TransformOperatorTest.
         */
         TransformOperatorTest(SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransformOperatorTest() {};

         /**
          * @brief Write grid points to file
          */
         void writeSetup();

         /**
          * @brief Write grid points to file
          */
         void writeGrid();

         /**
          * @brief Write weights to file
          */
         void writeWeights();

         /**
          * @brief Write projectors to file
          */
         virtual void writeProjectors() = 0;

         /**
          * @brief Write integrators to file
          */
         virtual void writeIntegrators() = 0;
         
      protected:
         /**
          * @brief Basename of the grid points file
          */
         const std::string mSetupBasename;

         /**
          * @brief Basename of the grid points file
          */
         const std::string mGridBasename;

         /**
          * @brief Basename of the weights file
          */
         const std::string mWeightsBasename;

         /**
          * @brief Extension for the projector files
          */
         const std::string mProjExt;

         /**
          * @brief Extension for the integrator files
          */
         const std::string mIntgExt;

         /**
          * @brief Specific test file prefix
          */
         std::string mPrefix;

         /**
          * @brief Ordering dependent file separation string
          */
         std::string mOrderString;

         /**
          * @brief Grid dimension
          */
         int mGridSize;

         /**
          * @brief Spectrum dimension
          */
         int mSpectrumSize;

         /**
          * @brief Maximum ordering index
          */
         int mOrderingN;

         /**
          * @brief Operator basis object
          */
         TBasis   mBasis;

         /**
          * @brief Get spectrum size depending on ordering index
          */
         virtual int spectrumSize(const int index) = 0;

         /**
          * @brief Write operator data to file
          *
          * @param name Name of the operator
          * @param ext File extension to use
          */
         template<class T, typename TExp> void writeOperator(std::string name, const PolynomialOperator<TExp>& (T::*Tptr)() const, std::string ext);


      private:
   };

   template <typename TBasis> TransformOperatorTest<TBasis>::TransformOperatorTest(SmartTruncation pTrunc)
      : mSetupBasename("setup.dat"), mGridBasename("grid.data"), mWeightsBasename("weights.data"), mProjExt(".proj"), mIntgExt(".intg"), mPrefix(""), mOrderString(""), mGridSize(0), mSpectrumSize(0), mOrderingN(0), mBasis(pTrunc)
   {
   }

   template <typename TBasis> void TransformOperatorTest<TBasis>::writeSetup()
   {
      // Create file object
      std::ofstream  file;

      // Filename base string
      std::string filename;

      // Create base name
      filename = this->mPrefix + this->mSetupBasename;

      // Open file
      file.open(filename.c_str());

      // write max order index to file
      file << this->mOrderingN << std::endl;

      // write Radial truncation information to file
      file << this->mGridSize << std::endl;
      file << this->mSpectrumSize;

      // close file
      file.close();
   }

   template <typename TBasis> void TransformOperatorTest<TBasis>::writeGrid()
   {
      // Create file object
      std::ofstream  file;

      // Filename base string
      std::string basename;

      // Stringstream for the conversion
      std::stringstream converter;

      // Convert grid size
      converter << this->mGridSize;

      // Create base name
      basename = this->mPrefix + "N" + converter.str() + "_" + this->mGridBasename;

      // Empty string stream
      converter.str("");

      // Open file
      file.open(basename.c_str(), std::ios::binary);

      // Write data to file
      file.write((char *)this->mBasis.at(0).grid().data(), this->mGridSize*sizeof(EPMFloat));

      // close file
      file.close();
   }

   template <typename TBasis> void TransformOperatorTest<TBasis>::writeWeights()
   {
      // Create file object
      std::ofstream  file;

      // Filename base string
      std::string basename;

      // Stringstream for the conversion
      std::stringstream converter;

      // Convert grid size
      converter << this->mGridSize;

      // Create base name
      basename = this->mPrefix + "N" + converter.str() + "_" + this->mWeightsBasename;

      // Empty string stream
      converter.str("");

      // Open file
      file.open(basename.c_str(), std::ios::binary);

      // Write data to file
      file.write((char *)this->mBasis.at(0).weights().data(), this->mGridSize*sizeof(EPMFloat));

      // close file
      file.close();
   }

   template <typename TBasis> template<class T, typename TExp> void TransformOperatorTest<TBasis>::writeOperator(std::string name, const PolynomialOperator<TExp>& (T::*Tptr)() const, std::string ext)
   {
      // Create file object
      std::ofstream  file;

      // Filename base string
      std::string basename;

      // Filename string
      std::string filename;

      // Stringstream for the conversion
      std::stringstream converter;

      // Convert grid size
      converter << this->mGridSize;

      // Create base name
      basename = this->mPrefix + "N" + converter.str() + "_" + name + this->mOrderString;

      // Empty string stream
      converter.str("");

      // Loop over all harmonic orders
      for(int i = 0; i < this->mOrderingN; ++i)
      {
         // Convert harmonic order
         converter << i;

         // Build filename
         filename = basename + converter.str() + ext;

         // Open file
         file.open(filename.c_str(), std::ios::binary);

         // Write data to file
         file.write((char*)(this->mBasis.at(i).*Tptr)().productOp().data(), this->mGridSize*this->spectrumSize(i)*sizeof(EPMFloat));

         // close file
         file.close();

         // empty stringstream
         converter.str("");
      }
   }

}

#endif // TRANSFORMOPERATORTEST_HPP
