/** \file DipoleFile.hpp
 *  \brief Implementation of a dipole diagnostic file
 */

#ifndef DIPOLEFILE_HPP
#define DIPOLEFILE_HPP

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
#include "IO/ASCII/DipoleFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a dipole diagnostics file
    *
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class DipoleFile: public ASCIIFieldWriter<TFieldType, ASCIIEWriter>
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
         DipoleFile(TFieldType& var, std::string name, const TimestepParameters &tsParams, const RadialBasisType &basis, const int format = 0);

         /**
         * @brief Destructor
         */
         virtual ~DipoleFile() {};

         /**
          * @brief Write energy to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;

         /**
          * @brief Output format flag
          */
         const int mFormat;

         /**
          * @brief Storage for the Dipole field
          */
         Array mDipole;
         /**
          * @brief Reference to the radial basis
          */
         const RadialBasisType&   mrBasis;

   };

   template <typename TFieldType> DipoleFile<TFieldType>::DipoleFile(TFieldType &var, std::string name, const TimestepParameters &tsParams, const RadialBasisType &basis, int format)
      : ASCIIFieldWriter<TFieldType, ASCIIEWriter>(var, name + DipoleFileDefs::BASENAME, DipoleFileDefs::EXTENSION, DipoleFileDefs::HEADER, DipoleFileDefs::VERSION), mrTSParams(tsParams), mFormat(format), mDipole(3), mrBasis(basis)
   {
   }

   template <typename TFieldType> void DipoleFile<TFieldType>::write()
   {
      // Compute dipole field
      //this->mrVar.rOc().updateSpectra();
      this->mDipole = this->mrVar.rOc().rPerturbation().dipole(this->mrBasis);

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write dipole
         this->mFile << std::setprecision(16) << this->mrTSParams.time() << "  " << this->mDipole.transpose();

         // Add newline at the end of ouput
         this->mFile << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

}

#endif // DIPOLEFILE_HPP
