/** \file PhysicalTorPolBase.hpp
 *  \brief Base of the implementation of Toroidal/Poloidal expanded field field
 */

#ifndef PHYSICALTORPOLBASE_HPP
#define PHYSICALTORPOLBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalFieldBase.hpp"
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
   template <typename TSimType> class PhysicalTorPolBase : public PhysicalFieldBase<TSimType>
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
         PhysicalTorPolBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalTorPolBase() {};

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
         
      protected:
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

   template<typename TSimType> PhysicalTorPolBase<TSimType>::PhysicalTorPolBase(SmartTruncation pTrunc, typename PhysicalTorPolBase<TSimType>::TransformType &transform)
      : PhysicalFieldBase<TSimType>(pTrunc, transform), mNeedCurlTransform(0), mPerturbation(pTrunc)
   {
   }

   template<typename TSimType> inline const TorPolField<TSimType>& PhysicalTorPolBase<TSimType>::perturbation() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline const TorPolField<TSimType>& PhysicalTorPolBase<TSimType>::totalField() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline TorPolField<TSimType>& PhysicalTorPolBase<TSimType>::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedCurlTransform = 0;

      return mPerturbation;
   }

}

#endif // PHYSICALTORPOLBASE_HPP
