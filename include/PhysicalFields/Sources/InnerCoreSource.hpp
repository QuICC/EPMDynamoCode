/** \file InnerCoreSource.hpp
 *  \brief Implementation of a field with a source term in inner core
 */

#ifndef INNERCORESOURCE_HPP
#define INNERCORESOURCE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a field with a source term in inner core
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSourceTraits Traits describing the field
    */
   template <typename TSimType, template <typename> class TSourceTraits> class InnerCoreSource: public TSourceTraits<TSimType>::DomainType
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from the spectral field from field traits
         typedef typename TSourceTraits<TSimType>::SpectralFieldType    SpectralFieldType;

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
         
      protected:

      private:
         /**
          * @brief Storage for Source term
          */
         SpectralFieldType mIcSrc;

   };

   template <typename TSimType, template <typename>  class TSourceTraits> InnerCoreSource<TSimType, TSourceTraits>::InnerCoreSource(SmartTruncation pTrunc, typename InnerCoreSource<TSimType,TSourceTraits>::TransformType &transform)
      : TSourceTraits<TSimType>::DomainType(pTrunc, transform), mIcSrc(pTrunc, true)
   {
   }

   template <typename TSimType, template <typename>  class TSourceTraits> const typename InnerCoreSource<TSimType, TSourceTraits>::SpectralFieldType& InnerCoreSource<TSimType, TSourceTraits>::icSrc() const
   {
      return this->mIcSrc;
   }

   template <typename TSimType, template <typename>  class TSourceTraits> typename InnerCoreSource<TSimType, TSourceTraits>::SpectralFieldType& InnerCoreSource<TSimType, TSourceTraits>::rIcSrc()
   {
      return this->mIcSrc;
   }

}

#endif // INNERCORESOURCE_HPP
