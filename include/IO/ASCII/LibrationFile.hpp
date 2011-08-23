/** \file LibrationFile.hpp
 *  \brief Implementation of a libration diagnostic file
 */

#ifndef LIBRATIONFILE_HPP
#define LIBRATIONFILE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IO/ASCII/LibrationFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a libration diagnostics file
    *
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class LibrationFile: public ASCIIFieldWriter<TFieldType, ASCIIEWriter>
   {
      public:
         /// Typedef for the Spectral radial Transform data type
         typedef typename SimulationConfig::NumericalScheme::RadialTransformType::PolynomialBasis RadialBasisType;

         /**
         * @brief Constructor
         *
         * @param var Field variable
         * @param name File name
         * @param tsParams Timestep parameters
         */
         LibrationFile(TFieldType& var, std::string name, const TimestepParameters &tsParams, const RadialBasisType &basis, const int format = 0);

         /**
         * @brief Destructor
         */
         virtual ~LibrationFile() {};

         /**
          * @brief Write energy to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Compute the poincare mode amplitude
          */
         void computePoincareMode();

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;

         /**
          * @brief Output format flag
          */
         const int mFormat;

         /**
          * @brief Storage for the Poincare mode amplitude
          */
         Array mAmplitude;

         /**
          * @brief Reference to the radial basis
          */
         const RadialBasisType&   mrBasis;

   };

   template <typename TFieldType> LibrationFile<TFieldType>::LibrationFile(TFieldType &var, std::string name, const TimestepParameters &tsParams, const RadialBasisType &basis, int format)
      : ASCIIFieldWriter<TFieldType, ASCIIEWriter>(var, name + LibrationFileDefs::BASENAME, LibrationFileDefs::EXTENSION, LibrationFileDefs::HEADER, LibrationFileDefs::VERSION), mrTSParams(tsParams), mFormat(format), mAmplitude(3), mrBasis(basis)
   {
   }

   template <typename TFieldType> void LibrationFile<TFieldType>::write()
   {
      // Compute poincare mode
      this->computePoincareMode();

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write energy
         this->mFile << this->mrTSParams.time() << "  " << this->mAmplitude.transpose();

         // Add newline at the end of ouput
         this->mFile << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

   template <typename TFieldType> void LibrationFile<TFieldType>::computePoincareMode()
   {
      this->mAmplitude = this->mrVar.rOc().rPerturbation().computeXYSolidComputation(this->mrBasis);
   }

}

#endif // LIBRATIONFILE_HPP
