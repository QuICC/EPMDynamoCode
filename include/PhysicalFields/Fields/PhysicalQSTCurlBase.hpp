/** \file PhysicalQSTCurlBase.hpp
 *  \brief Implementation of the base of the QST expanded physical field
 */

#ifndef PHYSICALQSTCURLBASE_HPP
#define PHYSICALQSTCURLBASE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "PhysicalFields/Fields/PhysicalRTPFieldCurl.hpp"
#include "GeneralFields/QSTField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the base of the QST expanded physical field
    *
    * The used spectral expansion is an QST expansion as it is a completely general
    * expansion
    */
   class PhysicalQSTCurlBase : public PhysicalRTPFieldCurl
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef QSTField  SpectralFieldType;

         /**
         * @brief Constructs underlying rtp and spectral fields
         *
         * @param pTrunc Truncation information
         * @param transform Transform object
         */
         PhysicalQSTCurlBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PhysicalQSTCurlBase() {};

         /**
          * @brief Get QST decomposition of the perturbation part
          */
         const QSTField&  perturbation() const;

         /**
          * @brief Get QST decomposition of the total field
          */
         const QSTField&  totalField() const;

         /**
          * @brief Set QST decomposition of the perturbation part
          */
         QSTField&  rPerturbation();

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();

      protected:
         /**
          * @brief Flag to check if curl transform has already been
          *          computed
          */
         int  mNeedCurlTransform;

         /**
          * @brief Spectral QST decomposition of the field
          */
         QSTField    mPerturbation;

      private:
   };

   inline const QSTField& PhysicalQSTCurlBase::perturbation() const
   {
      return this->mPerturbation;
   }

   inline const QSTField& PhysicalQSTCurlBase::totalField() const
   {
      return this->mPerturbation;
   }

   inline QSTField& PhysicalQSTCurlBase::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedCurlTransform = 0;
      
      return this->mPerturbation;
   }

}

#endif // PHYSICALQSTCURLBASE_HPP
