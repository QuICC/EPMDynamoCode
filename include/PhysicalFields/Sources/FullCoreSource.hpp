/** \file FullCoreSource.hpp
 *  \brief Implementation of a field with a source term in full core
 */

#ifndef FULLCORESOURCE_HPP
#define FULLCORESOURCE_HPP

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

namespace EPMDynamo {

   /**
    * \brief Implementation of a field with a source term in full core
    *
    * \tparam TSourceTraits Traits describing the field
    */
   template <typename TSourceTraits> class FullCoreSource: public TSourceTraits::DomainType
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef from the spectral field from field traits
         typedef typename TSourceTraits::SpectralFieldType    SpectralFieldType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Transform object
          */
         FullCoreSource(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~FullCoreSource() {};

         /**
          * @brief Get values of inner core source term
          */
         const SpectralFieldType& icSrc() const;

         /**
          * @brief Set values of inner core source term
          */
         SpectralFieldType&  rIcSrc();

         /**
          * @brief Get values of outer core source term
          */
         const SpectralFieldType& ocSrc() const;

         /**
          * @brief Set values of outer core source term
          */
         SpectralFieldType&  rOcSrc();

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:

      private:
         /**
          * @brief Storage for inner core source term
          */
         SpectralFieldType mIcSrc;

         /**
          * @brief Storage for outer core source term
          */
         SpectralFieldType mOcSrc;

   };

   template <typename TSourceTraits> FullCoreSource<TSourceTraits>::FullCoreSource(SmartTruncation pTrunc, typename FullCoreSource<TSourceTraits>::TransformType &transform)
      : TSourceTraits::DomainType(pTrunc, transform), mIcSrc(pTrunc, true), mOcSrc(pTrunc, true)
   {
   }

   template <typename TSourceTraits> const typename FullCoreSource<TSourceTraits>::SpectralFieldType& FullCoreSource<TSourceTraits>::icSrc() const
   {
      return this->mIcSrc;
   }

   template <typename TSourceTraits> typename FullCoreSource<TSourceTraits>::SpectralFieldType& FullCoreSource<TSourceTraits>::rIcSrc()
   {
      return this->mIcSrc;
   }

   template <typename TSourceTraits> const typename FullCoreSource<TSourceTraits>::SpectralFieldType& FullCoreSource<TSourceTraits>::ocSrc() const
   {
      return this->mOcSrc;
   }

   template <typename TSourceTraits> typename FullCoreSource<TSourceTraits>::SpectralFieldType& FullCoreSource<TSourceTraits>::rOcSrc()
   {
      return this->mOcSrc;
   }

   template <typename TSourceTraits> void FullCoreSource<TSourceTraits>::initialiseZeros()
   {
      TSourceTraits::DomainType::initialiseZeros();

      this->mIcSrc.intialiseZeros();

      this->mOcSrc.intialiseZeros();
   }

}

#endif // FULLCORESOURCE_HPP
