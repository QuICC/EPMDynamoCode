/** \file PhysicalScalarBase.hpp
 *  \brief Implementation of physical scalar base
 */

#ifndef PHYSICALSCALARBASE_HPP
#define PHYSICALSCALARBASE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPScalar.hpp"
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of physical scalar base
    */
   class PhysicalScalarBase: public PhysicalRTPScalar
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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:
         /**
          * @brief Need any transform?
          */
         bool needAnyTransform() const;

         /**
          * @brief Spectral expansion scalar
          */
         ScalarType   mPerturbation;

      private:
   };

   inline bool PhysicalScalarBase::needAnyTransform() const
   {
      return (this->mNeedTransform == 0);
   }

   inline const typename PhysicalScalarBase::ScalarType& PhysicalScalarBase::perturbation() const
   {
      return this->mPerturbation;
   }

   inline const typename PhysicalScalarBase::ScalarType& PhysicalScalarBase::totalField() const
   {
      return this->mPerturbation;
   }

   inline typename PhysicalScalarBase::ScalarType& PhysicalScalarBase::rPerturbation()
   {
      this->mNeedTransform = 0;

      return this->mPerturbation;
   }

}

#endif // PHYSICALSCALARBASE_HPP
