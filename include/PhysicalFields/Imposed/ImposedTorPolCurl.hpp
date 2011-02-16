/** \file ImposedTorPolCurl.hpp
 *  \brief Implementation of Toroidal/Poloidal expanded imposed field
 */

#ifndef IMPOSEDTORPOLCURL_HPP
#define IMPOSEDTORPOLCURL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolCurlBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of Toroidal/Poloidal expanded imposed field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ImposedTorPolCurl : public PhysicalTorPolCurlBase<TSimType>
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

   template<typename TSimType> ImposedTorPolCurl<TSimType>::ImposedTorPolCurl(SmartTruncation pTrunc, typename ImposedTorPolCurl<TSimType>::TransformType &transform)
      : PhysicalTorPolCurlBase<TSimType>(pTrunc, transform), mTotalField(pTrunc), mImposedField(pTrunc)
   {
   }

   template<typename TSimType> inline const typename ImposedTorPolCurl<TSimType>::SpectralFieldType& ImposedTorPolCurl<TSimType>::totalField() const
   {
      return this->mTotalField;
   }

   template<typename TSimType> inline const typename ImposedTorPolCurl<TSimType>::SpectralFieldType& ImposedTorPolCurl<TSimType>::imposed() const
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline typename ImposedTorPolCurl<TSimType>::SpectralFieldType& ImposedTorPolCurl<TSimType>::rImposed()
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline void ImposedTorPolCurl<TSimType>::updateTotalField()
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

   template<typename TSimType> void ImposedTorPolCurl<TSimType>::initialiseZeros()
   {
      PhysicalTorPolCurlBase<TSimType>::initialiseZeros();

      this->mImposedField.initialiseZeros();
   }

}

#endif // IMPOSEDTORPOLCURL_HPP
