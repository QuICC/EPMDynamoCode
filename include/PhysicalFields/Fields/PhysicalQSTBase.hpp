/** \file PhysicalQSTBase.hpp
 *  \brief Implementation of the base of the QST expanded physical field
 */

#ifndef PHYSICALQSTBASE_HPP
#define PHYSICALQSTBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"
#include "PhysicalFields/Fields/PhysicalFieldBase.hpp"
#include "GeneralFields/QSTField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the base of the QST expanded physical field
    *
    * The used spectral expansion is an QST expansion as it is a completely general
    * expansion
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalQSTBase : public PhysicalFieldBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef QSTField<TSimType>  SpectralFieldType;

         /**
         * @brief Constructs underlying rtp and spectral fields
         *
         * @param pTrunc Truncation information
         * @param transform Transform object
         */
         PhysicalQSTBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PhysicalQSTBase() {};

         /**
          * @brief Get QST decomposition of the perturbation part
          */
         const QSTField<TSimType>&  perturbation() const;

         /**
          * @brief Get QST decomposition of the total field
          */
         const QSTField<TSimType>&  totalField() const;

         /**
          * @brief Set QST decomposition of the perturbation part
          */
         QSTField<TSimType>&  rPerturbation();

      protected:
         /**
          * @brief Flag to check if curl transform has already been
          *          computed
          */
         int  mNeedCurlTransform;

         /**
          * @brief Spectral QST decomposition of the field
          */
         QSTField<TSimType>    mPerturbation;

      private:
   };

   template<typename TSimType> PhysicalQSTBase<TSimType>::PhysicalQSTBase(SmartTruncation pTrunc, typename PhysicalQSTBase<TSimType>::TransformType &transform)
      : PhysicalFieldBase<TSimType>(pTrunc, transform), mNeedCurlTransform(0), mPerturbation(pTrunc)
   {
   }

   template<typename TSimType> inline const QSTField<TSimType>& PhysicalQSTBase<TSimType>::perturbation() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline const QSTField<TSimType>& PhysicalQSTBase<TSimType>::totalField() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline QSTField<TSimType>& PhysicalQSTBase<TSimType>::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedCurlTransform = 0;
      
      return this->mPerturbation;
   }

}

#endif // PHYSICALQSTBASE_HPP
