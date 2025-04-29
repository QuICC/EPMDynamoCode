/** \file AngularMomentumFile.hpp
 *  \brief Implementation of a angular momentum diagnostic file
 */

#ifndef ANGULARMOMENTUMFILE_HPP
#define ANGULARMOMENTUMFILE_HPP

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
#include "IO/ASCII/AngularMomentumFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a angular momentum diagnostics file
    *
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class AngularMomentumFile: public ASCIIFieldWriter<TFieldType, ASCIIEWriter>
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
         AngularMomentumFile(TFieldType& var, std::string name, const TimestepParameters &tsParams, const RadialBasisType &basis, const int format = 0);

         /**
         * @brief Destructor
         */
         virtual ~AngularMomentumFile() {};

         /**
          * @brief Write angular momentum to file
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
          * @brief Storage for the Poincare mode amplitude
          */
         Array mMomentum;

         /**
          * @brief Reference to the radial basis
          */
         const RadialBasisType&   mrBasis;

   };

   template <typename TFieldType> AngularMomentumFile<TFieldType>::AngularMomentumFile(TFieldType &var, std::string name, const TimestepParameters &tsParams, const RadialBasisType &basis, int format)
      : ASCIIFieldWriter<TFieldType, ASCIIEWriter>(var, name + AngularMomentumFileDefs::BASENAME, AngularMomentumFileDefs::EXTENSION, AngularMomentumFileDefs::HEADER, AngularMomentumFileDefs::VERSION), mrTSParams(tsParams), mFormat(format), mMomentum(4), mrBasis(basis)
   {
   }

   template <typename TFieldType> void AngularMomentumFile<TFieldType>::write()
   {
      Array corr = Array::Zero(3);
      this->mMomentum = this->mrVar.rOc().rPerturbation().computeXYZAngularMomentum(this->mrBasis);

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write angular momentum
         this->mFile << std::setprecision(16) << this->mrTSParams.time() << "  " << this->mMomentum.transpose();

         // Add newline at the end of ouput
         this->mFile << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

}

#endif // ANGULARMOMENTUMFILE_HPP
