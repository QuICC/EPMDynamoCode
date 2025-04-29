/** \file ImposedScalar.hpp
 *  \brief Implementation of scalar imposed field
 */

#ifndef IMPOSEDSCALAR_HPP
#define IMPOSEDSCALAR_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalScalarBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of scalar imposed field
    */
   class ImposedScalar : public PhysicalScalarBase
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef PhysicalScalarBase::SpectralFieldType  SpectralFieldType;

         /**
          * @brief Constructs the underlying rtp and spectral fields
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         ImposedScalar(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~ImposedScalar() {};

         /**
          * @brief Get decomposition of the field (total field)
          */
         const SpectralFieldType&  totalField() const;

         /**
          * @brief Get imposed field
          */
         const SpectralFieldType&  imposed() const;

         /**
          * @brief set imposed field 
          */
         SpectralFieldType&  rImposed();

         /**
          * @brief Update the value of the total field (perturbation + imposed)
          */
         virtual void updateTotalField();

         /**
          * @brief Initialise to zeros
          */
         void initialiseZeros();
         
      protected:
         /**
          * @brief Spectral decomposition of the field (total field)
          */
         SpectralFieldType    mTotalField;

         /**
          * @brief Spectral decomposition of the imposed field
          */
         SpectralFieldType    mImposedField;

      private:
   };

   inline const ImposedScalar::SpectralFieldType& ImposedScalar::totalField() const
   {
      return this->mTotalField;
   }

   inline const ImposedScalar::SpectralFieldType& ImposedScalar::imposed() const
   {
      return this->mImposedField;
   }

   inline ImposedScalar::SpectralFieldType& ImposedScalar::rImposed()
   {
      return this->mImposedField;
   }

   inline void ImposedScalar::updateTotalField()
   {
      if(this->needAnyTransform())
      {
         // Set total scalar
         int nL = this->perturbation().nL();
         for(int l = this->perturbation().minL(); l < nL ; ++l)
         {
            this->mTotalField.rLShell(l) = this->perturbation().lshell(l) + this->imposed().lshell(l);
         }
      }
   }

}

#endif // IMPOSEDSCALAR_HPP
