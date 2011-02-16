/** \file ImposedScalar.hpp
 *  \brief Implementation of scalar imposed field
 */

#ifndef IMPOSEDSCALAR_HPP
#define IMPOSEDSCALAR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalScalarBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of scalar imposed field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ImposedScalar : public PhysicalScalarBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef typename PhysicalScalarBase<TSimType>::SpectralFieldType  SpectralFieldType;

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

   template<typename TSimType> ImposedScalar<TSimType>::ImposedScalar(SmartTruncation pTrunc, typename ImposedScalar<TSimType>::TransformType &transform)
      : PhysicalScalarBase<TSimType>(pTrunc, transform), mTotalField(pTrunc, true), mImposedField(pTrunc, true)
   {
   }

   template<typename TSimType> inline const typename ImposedScalar<TSimType>::SpectralFieldType& ImposedScalar<TSimType>::totalField() const
   {
      return this->mTotalField;
   }

   template<typename TSimType> inline const typename ImposedScalar<TSimType>::SpectralFieldType& ImposedScalar<TSimType>::imposed() const
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline typename ImposedScalar<TSimType>::SpectralFieldType& ImposedScalar<TSimType>::rImposed()
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline void ImposedScalar<TSimType>::updateTotalField()
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

   template<typename TSimType> void ImposedScalar<TSimType>::initialiseZeros()
   {
      PhysicalScalarBase<TSimType>::initialiseZeros();

      this->mImposedField.initialiseZeros();
   }

}

#endif // IMPOSEDSCALAR_HPP
