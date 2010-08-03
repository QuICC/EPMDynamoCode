/** \file PhysicalScalarBase.hpp
 *  \brief Implementation of physical scalar base
 */

#ifndef PHYSICALSCALARBASE_HPP
#define PHYSICALSCALARBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalScalarFieldBase.hpp"
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of physical scalar base
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalScalarBase: public PhysicalScalarFieldBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType  ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType  TransformType;

         /// Typedef for the spectral field type
         typedef ScalarType SpectralFieldType;

         /**
         * @brief Constructor
         *
         * @param pTrunc Truncation information
         * @param transform Reference to transform object
         */
         PhysicalScalarBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalScalarBase() {};

         /**
          * @brief Get spectral expansion scalar (perturbation part)
          */
         const ScalarType&  perturbation() const;

         /**
          * @brief Get spectral expansion scalar (Total scalar field)
          */
         const ScalarType&  totalField() const;

         /**
          * @brief Set spectral expansion scalar
          */
         ScalarType&  rPerturbation();
         
      protected:
         /**
          * @brief Flag to check if grad transform has already been computed
          */
         int  mNeedGradTransform;

         /**
          * @brief Spectral expansion scalar
          */
         ScalarType   mPerturbation;

      private:
   };
   
   template <typename TSimType> PhysicalScalarBase<TSimType>::PhysicalScalarBase(SmartTruncation pTrunc, typename PhysicalScalarBase<TSimType>::TransformType &transform)
      : PhysicalScalarFieldBase<TSimType>(pTrunc, transform), mNeedGradTransform(0), mPerturbation(pTrunc, true)
   {
   }

   template <typename TSimType> inline const typename PhysicalScalarBase<TSimType>::ScalarType& PhysicalScalarBase<TSimType>::perturbation() const
   {
      return this->mPerturbation;
   }

   template <typename TSimType> inline const typename PhysicalScalarBase<TSimType>::ScalarType& PhysicalScalarBase<TSimType>::totalField() const
   {
      return this->mPerturbation;
   }

   template <typename TSimType> inline typename PhysicalScalarBase<TSimType>::ScalarType& PhysicalScalarBase<TSimType>::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedGradTransform = 0;

      return this->mPerturbation;
   }

}

#endif // PHYSICALSCALARBASE_HPP
