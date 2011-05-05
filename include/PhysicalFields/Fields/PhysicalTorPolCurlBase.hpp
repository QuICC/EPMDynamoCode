/** \file PhysicalTorPolCurlBase.hpp
 *  \brief Base of the implementation of Toroidal/Poloidal expanded field field
 *  with curl
 */

#ifndef PHYSICALTORPOLCURlBASE_HPP
#define PHYSICALTORPOLCURlBASE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPFieldCurl.hpp"
#include "GeneralFields/TorPolField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of Toroidal/Poloidal expanded field field
    *
    * The used spectral expansion is a Toroidal/Poloidal expansion to take advantage of
    * the divergence free aspect of the field
    */
   class PhysicalTorPolCurlBase : public PhysicalRTPFieldCurl
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
         PhysicalTorPolCurlBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalTorPolCurlBase() {};

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
          * @brief Flag to check if curl transform has already been computed
          */
         int  mNeedCurlTransform;

         /**
          * @brief Spectral toroidal/poloidal decomposition of the field
          */
         TorPolField    mPerturbation;

      private:
   };

   inline bool PhysicalTorPolCurlBase::needAnyTransform() const
   {
      return ((this->mNeedTransform == 0) && (this->mNeedCurlTransform == 0));
   }

   inline const TorPolField& PhysicalTorPolCurlBase::perturbation() const
   {
      return this->mPerturbation;
   }

   inline const TorPolField& PhysicalTorPolCurlBase::totalField() const
   {
      return this->mPerturbation;
   }

   inline TorPolField& PhysicalTorPolCurlBase::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedCurlTransform = 0;

      return this->mPerturbation;
   }

}

#endif // PHYSICALTORPOLCURlBASE_HPP
