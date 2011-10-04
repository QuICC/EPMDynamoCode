/** \file TorPolField.hpp
 *  \brief Implementation of the Toroidal/Poloidal vector field decomposition
 */

#ifndef TORPOLFIELD_HPP
#define TORPOLFIELD_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/MathConstants.hpp"
#include "GeneralFields/SpectralFieldBase.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Toroidal/Poloidal vector field decomposition
    */
   class TorPolField: public SpectralFieldBase
   {
      /// Typedef from Simulation trait to local truncation type
      typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

      /// Typedef from Simulation trait to local truncation type
      typedef SimulationConfig::NumericalScheme::RadialBasisType    RadialBasisType;

      public:
         /**
         * @brief Constructs the field using the given spectral truncation
         *
         * @param pTrunc Smart pointer to truncation information
         * @param needsL0  Needs to consider L=0 mode ?
         */
         TorPolField(SmartTruncation pTrunc, const bool needsL0 = false);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TorPolField() {};

         /**
          * @brief Get the Toroidal component of the field
          */
         const ScalarType&  tor() const;

         /**
          * @brief Set the Toroidal component of the field
          */
         ScalarType&  rTor();

         /**
          * @brief Get the Poloidal component of the field
          */
         const ScalarType&  pol() const;

         /**
          * @brief Set the Poloidal component of the field
          */
         ScalarType&  rPol();

         /**
          * @brief Compute the toroidal energy spectra
          *
          * @param radBasis Radial basis 
          */
         void computeTorSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Compute the poloidal energy spectra
          *
          * @param radBasis Radial basis 
          */
         void computePolSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Compute the \f$\hat{x}\f$ and \f$\hat{y}\f$ solid body projections
          *
          * @param radBasis Radial basis 
          */
         Array computeXYSolidProjection(const RadialBasisType &radBasis);

         /**
          * @brief Get the energy in the toroidal component
          */
         Array  torEnergy() const;

         /**
          * @brief Get the energy in the poloidal component
          */
         Array  polEnergy() const;

         /**
          * @brief Get the energy spectrum in the toroidal component
          */
         const Array&  torSpectrumL() const;

         /**
          * @brief Get the energy spectrum in the poloidal component
          */
         const Array&  polSpectrumL() const;

         /**
          * @brief Get the energy spectrum in the toroidal component
          */
         const Array&  torSpectrumM() const;

         /**
          * @brief Get the energy spectrum in the poloidal component
          */
         const Array&  polSpectrumM() const;

         /**
          * @brief Set the normalisation factors
          *
          * @param eFactor Energy normalisation factor
          */
         void setEnergyScale(const EPMFloat&  eFactor);

         /**
          * @brief Rescale the Toroidal/Poloidal components by a constant
          *
          * @param scale Scaling factor
          */
         void rescale(const EPMFloat scale);

         /**
          * @brief initialise field to zeros
          */
         void initialiseZeros();

         /**
          * @brief Check field set during initialisation
          */
         void checkInitialisation();

      protected:

      private:
         /**
          * @brief Storage for the Toroidal component
          */
         ScalarType   mTor;

         /**
          * @brief Storage for the Poloidal component
          */
         ScalarType   mPol;
   };

   inline const TorPolField::ScalarType& TorPolField::tor() const
   {
      return this->mTor;
   }

   inline TorPolField::ScalarType& TorPolField::rTor()
   {
      return this->mTor;
   }

   inline const TorPolField::ScalarType& TorPolField::pol() const
   {
      return this->mPol;
   }

   inline TorPolField::ScalarType& TorPolField::rPol()
   {
      return this->mPol;
   }
}

#endif // TORPOLFIELD_HPP
