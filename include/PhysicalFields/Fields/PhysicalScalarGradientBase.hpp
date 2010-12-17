/** \file PhysicalScalarGradientBase.hpp
 *  \brief Implementation of physical scalar base with gradient
 */

#ifndef PHYSICALSCALARGRADIENTBASE_HPP
#define PHYSICALSCALARGRADIENTBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalRTPScalarGradient.hpp"
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of physical scalar base
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalScalarGradientBase: public PhysicalRTPScalarGradient<TSimType>
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
         PhysicalScalarGradientBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalScalarGradientBase() {};

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
          * @brief Need any transform ?
          *
          * This is important to generalise cases with or without gradient
          */
         bool needAnyTransform() const;

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
   
   template <typename TSimType> PhysicalScalarGradientBase<TSimType>::PhysicalScalarGradientBase(SmartTruncation pTrunc, typename PhysicalScalarGradientBase<TSimType>::TransformType &transform)
      : PhysicalRTPScalarGradient<TSimType>(pTrunc, transform), mNeedGradTransform(0), mPerturbation(pTrunc, true)
   {
   }

   template <typename TSimType> inline bool PhysicalScalarGradientBase<TSimType>::needAnyTransform() const
   {
      return ((this->mNeedTransform == 0) && (this->mNeedGradTransform == 0));
   }

   template <typename TSimType> inline const typename PhysicalScalarGradientBase<TSimType>::ScalarType& PhysicalScalarGradientBase<TSimType>::perturbation() const
   {
      return this->mPerturbation;
   }

   template <typename TSimType> inline const typename PhysicalScalarGradientBase<TSimType>::ScalarType& PhysicalScalarGradientBase<TSimType>::totalField() const
   {
      return this->mPerturbation;
   }

   template <typename TSimType> inline typename PhysicalScalarGradientBase<TSimType>::ScalarType& PhysicalScalarGradientBase<TSimType>::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedGradTransform = 0;

      return this->mPerturbation;
   }

}

#endif // PHYSICALSCALARGRADIENTBASE_HPP
