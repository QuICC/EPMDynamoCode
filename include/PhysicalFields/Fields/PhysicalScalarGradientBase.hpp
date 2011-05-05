/** \file PhysicalScalarGradientBase.hpp
 *  \brief Implementation of physical scalar base with gradient
 */

#ifndef PHYSICALSCALARGRADIENTBASE_HPP
#define PHYSICALSCALARGRADIENTBASE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPScalarGradient.hpp"
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of physical scalar base
    */
   class PhysicalScalarGradientBase: public PhysicalRTPScalarGradient
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType  ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType  TransformType;

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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
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

   inline bool PhysicalScalarGradientBase::needAnyTransform() const
   {
      return ((this->mNeedTransform == 0) && (this->mNeedGradTransform == 0));
   }

   inline const typename PhysicalScalarGradientBase::ScalarType& PhysicalScalarGradientBase::perturbation() const
   {
      return this->mPerturbation;
   }

   inline const typename PhysicalScalarGradientBase::ScalarType& PhysicalScalarGradientBase::totalField() const
   {
      return this->mPerturbation;
   }

   inline typename PhysicalScalarGradientBase::ScalarType& PhysicalScalarGradientBase::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedGradTransform = 0;

      return this->mPerturbation;
   }

}

#endif // PHYSICALSCALARGRADIENTBASE_HPP
