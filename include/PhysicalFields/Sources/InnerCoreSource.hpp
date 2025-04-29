/** \file InnerCoreSource.hpp
 *  \brief Implementation of a field with a source term in inner core
 */

#ifndef INNERCORESOURCE_HPP
#define INNERCORESOURCE_HPP

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
    * \brief Implementation of a field with a source term in inner core
    *
    * \tparam TSourceTraits Traits describing the field
    */
   template <typename TSourceTraits> class InnerCoreSource: public TSourceTraits::DomainType
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
         InnerCoreSource(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~InnerCoreSource() {};

         /**
          * @brief Get values of inner core source term
          */
         const SpectralFieldType& icSrc() const;

         /**
          * @brief Set values of inner core source term
          */
         SpectralFieldType&  rIcSrc();

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:

      private:
         /**
          * @brief Storage for Source term
          */
         SpectralFieldType mIcSrc;

   };

   template <typename TSourceTraits> InnerCoreSource<TSourceTraits>::InnerCoreSource(SmartTruncation pTrunc, typename InnerCoreSource<TSourceTraits>::TransformType &transform)
      : TSourceTraits::DomainType(pTrunc, transform), mIcSrc(pTrunc, true)
   {
   }

   template <typename TSourceTraits> const typename InnerCoreSource<TSourceTraits>::SpectralFieldType& InnerCoreSource<TSourceTraits>::icSrc() const
   {
      return this->mIcSrc;
   }

   template <typename TSourceTraits> typename InnerCoreSource<TSourceTraits>::SpectralFieldType& InnerCoreSource<TSourceTraits>::rIcSrc()
   {
      return this->mIcSrc;
   }

   template <typename TSourceTraits>  void InnerCoreSource<TSourceTraits>::initialiseZeros()
   {
      TSourceTraits::DomainType::initialiseZeros();

      this->mIcSrc.initialiseZeros();
   }

}

#endif // INNERCORESOURCE_HPP
