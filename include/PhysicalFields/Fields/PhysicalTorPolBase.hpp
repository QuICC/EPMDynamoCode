/** \file PhysicalTorPolBase.hpp
 *  \brief Base of the implementation of Toroidal/Poloidal expanded field field
 */

#ifndef PHYSICALTORPOLBASE_HPP
#define PHYSICALTORPOLBASE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPField.hpp"
#include "GeneralFields/TorPolField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of Toroidal/Poloidal expanded field field
    *
    * The used spectral expansion is a Toroidal/Poloidal expansion to take advantage of
    * the divergence free aspect of the field
    */
   class PhysicalTorPolBase : public PhysicalRTPField
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef TorPolField  SpectralFieldType;

         /**
          * @brief Constructs the underlying rtp and spectral fields
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         PhysicalTorPolBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalTorPolBase() {};

         /**
          * @brief Get Toroidal/Poloidal decomposition of the field (perturbation part)
          */
         const TorPolField&  perturbation() const;

         /**
          * @brief Get Toroidal/Poloidal decomposition of the field (total field)
          */
         const TorPolField&  totalField() const;

         /**
          * @brief Set Toroidal/Poloidal decomposition of the field (perturbation part)
          */
         TorPolField&  rPerturbation();

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
          * @brief Spectral toroidal/poloidal decomposition of the field
          */
         TorPolField    mPerturbation;

      private:
   };

   inline bool PhysicalTorPolBase::needAnyTransform() const
   {
      return (this->mNeedTransform == 0);
   }

   inline const TorPolField& PhysicalTorPolBase::perturbation() const
   {
      return this->mPerturbation;
   }

   inline const TorPolField& PhysicalTorPolBase::totalField() const
   {
      return this->mPerturbation;
   }

   inline TorPolField& PhysicalTorPolBase::rPerturbation()
   {
      this->mNeedTransform = 0;

      return this->mPerturbation;
   }

}

#endif // PHYSICALTORPOLBASE_HPP
