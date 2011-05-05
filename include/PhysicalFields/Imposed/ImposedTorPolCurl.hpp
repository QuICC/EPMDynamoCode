/** \file ImposedTorPolCurl.hpp
 *  \brief Implementation of Toroidal/Poloidal expanded imposed field
 */

#ifndef IMPOSEDTORPOLCURL_HPP
#define IMPOSEDTORPOLCURL_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalTorPolCurlBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of Toroidal/Poloidal expanded imposed field
    */
   class ImposedTorPolCurl : public PhysicalTorPolCurlBase
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
         ImposedTorPolCurl(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~ImposedTorPolCurl() {};

         /**
          * @brief Get decomposition of the field (total field)
          */
         const SpectralFieldType&  totalField() const;

         /**
          * @brief Get imposed field
          */
         const SpectralFieldType&  imposed() const;

         /**
          * @brief set imposed field 
          */
         SpectralFieldType&  rImposed();

         /**
          * @brief Update the value of the total field (perturbation + imposed)
          */
         virtual void updateTotalField();

         /**
          * @brief Initialise to zeros
          */
         void initialiseZeros();
         
      protected:
         /**
          * @brief Spectral decomposition of the field (total field)
          */
         SpectralFieldType    mTotalField;

         /**
          * @brief Spectral decomposition of the imposed field
          */
         SpectralFieldType    mImposedField;

      private:
   };

   inline const ImposedTorPolCurl::SpectralFieldType& ImposedTorPolCurl::totalField() const
   {
      return this->mTotalField;
   }

   inline const ImposedTorPolCurl::SpectralFieldType& ImposedTorPolCurl::imposed() const
   {
      return this->mImposedField;
   }

   inline ImposedTorPolCurl::SpectralFieldType& ImposedTorPolCurl::rImposed()
   {
      return this->mImposedField;
   }

   inline void ImposedTorPolCurl::updateTotalField()
   {
      if(this->needAnyTransform())
      {
         // Set toroidal part
         int nL = this->perturbation().tor().nL();
         for(int l = this->perturbation().tor().minL(); l < nL ; ++l)
         {
            this->mTotalField.rTor().rLShell(l) = this->perturbation().tor().lshell(l) + this->imposed().tor().lshell(l);
         }

         // Set poloidal part
         nL = this->perturbation().pol().nL();
         for(int l = this->perturbation().pol().minL(); l < nL ; ++l)
         {
            this->mTotalField.rPol().rLShell(l) = this->perturbation().pol().lshell(l) + this->imposed().pol().lshell(l);
         }
      }
   }

}

#endif // IMPOSEDTORPOLCURL_HPP
