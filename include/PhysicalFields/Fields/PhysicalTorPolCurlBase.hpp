/** \file PhysicalTorPolCurlBase.hpp
 *  \brief Base of the implementation of Toroidal/Poloidal expanded field field
 *  with curl
 */

#ifndef PHYSICALTORPOLCURlBASE_HPP
#define PHYSICALTORPOLCURlBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalRTPFieldCurl.hpp"
#include "GeneralFields/TorPolField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of Toroidal/Poloidal expanded field field
    *
    * The used spectral expansion is a Toroidal/Poloidal expansion to take advantage of
    * the divergence free aspect of the field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalTorPolCurlBase : public PhysicalRTPFieldCurl<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef TorPolField<TSimType>  SpectralFieldType;

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
         const TorPolField<TSimType>&  perturbation() const;

         /**
          * @brief Get Toroidal/Poloidal decomposition of the field (total field)
          */
         const TorPolField<TSimType>&  totalField() const;

         /**
          * @brief Set Toroidal/Poloidal decomposition of the field (perturbation part)
          */
         TorPolField<TSimType>&  rPerturbation();

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
         TorPolField<TSimType>    mPerturbation;

      private:
   };

   template<typename TSimType> PhysicalTorPolCurlBase<TSimType>::PhysicalTorPolCurlBase(SmartTruncation pTrunc, typename PhysicalTorPolCurlBase<TSimType>::TransformType &transform)
      : PhysicalRTPFieldCurl<TSimType>(pTrunc, transform), mNeedCurlTransform(0), mPerturbation(pTrunc)
   {
   }

   template<typename TSimType> inline bool PhysicalTorPolCurlBase<TSimType>::needAnyTransform() const
   {
      return ((this->mNeedTransform == 0) && (this->mNeedCurlTransform == 0));
   }

   template<typename TSimType> inline const TorPolField<TSimType>& PhysicalTorPolCurlBase<TSimType>::perturbation() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline const TorPolField<TSimType>& PhysicalTorPolCurlBase<TSimType>::totalField() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline TorPolField<TSimType>& PhysicalTorPolCurlBase<TSimType>::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedCurlTransform = 0;

      return this->mPerturbation;
   }

   template<typename TSimType> void  PhysicalTorPolCurlBase<TSimType>::initialiseZeros()
   {
      PhysicalRTPFieldCurl<TSimType>::initialiseZeros();

      this->mPerturbation.initialiseZeros();
   }

}

#endif // PHYSICALTORPOLCURlBASE_HPP
