/** \file QSTField.hpp
 *  \brief Implementation of the QST field decomposition
 */

#ifndef QSTFIELD_HPP
#define QSTFIELD_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/SpectralFieldBase.hpp"
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the QST field decomposition
    */
   class QSTField: public SpectralFieldBase
   {
      /// Typedef from Simulation trait to local truncation type
      typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

      /// Typedef from Simulation trait to local truncation type
      typedef SimulationConfig::NumericalScheme::RadialBasisType    RadialBasisType;

      public:
         /**
         * @brief Constructs the field using the given spectral truncation information
         *
         * @param pTrunc Truncation information pointer
         */
         QSTField(SmartTruncation pTrunc);

         /**
         * @brief Simple empty constructor
         */
         virtual ~QSTField() {};

         /**
          * @brief Get Q component of the field
          */
         const ScalarType&  q() const;

         /**
          * @brief Set Q component of the field
          */
         ScalarType&  rQ();

         /**
          * @brief Get S component of the field
          */
         const ScalarType&  s() const;

         /**
          * @brief Set S component of the field
          */
         ScalarType&  rS();

         /**
          * @brief Get T component of the field
          */
         const ScalarType&  t() const;

         /**
          * @brief Set T component of the field
          */
         ScalarType&  rT();

         /**
          * @brief Rescale the QST components by a constant
          *
          * @param scale Scaling factor
          */
         void rescale(const EPMFloat scale);

         /**
          * @brief Define the energy non-dimensionlisation scales
          *
          * @param eFactor Energy multiplicative factor
          */
         void setEnergyScale(const EPMFloat&  eFactor);

         /**
          * @brief Compute the energy spectra of the Q QST Component
          *
          * \param radBasis Radial basis
          */
         void  computeQSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Compute the energy spectra of the S QST Component
          *
          * \param radBasis Radial basis
          */
         void  computeSSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Compute the energy spectra of the T QST Component
          *
          * \param radBasis Radial basis
          */
         void  computeTSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Initialise field to zero
          */
         void initialiseZeros();
         
      protected:

      private:
         /**
          * @brief Q component of the QST field
          */
         ScalarType   mQ;

         /**
          * @brief S component of the QST field
          */
         ScalarType   mS;

         /**
          * @brief T component of the QST field
          */
         ScalarType   mT;
   };

   inline const typename QSTField::ScalarType& QSTField::q() const
   {
      return this->mQ;
   }

   inline  typename QSTField::ScalarType& QSTField::rQ()
   {
      return this->mQ;
   }

   inline const  typename QSTField::ScalarType& QSTField::s() const
   {
      return this->mS;
   }

   inline  typename QSTField::ScalarType& QSTField::rS()
   {
      return this->mS;
   }

   inline const  typename QSTField::ScalarType& QSTField::t() const
   {
      return this->mT;
   }

   inline  typename QSTField::ScalarType& QSTField::rT()
   {
      return this->mT;
   }
}

#endif // QSTFIELD_HPP
