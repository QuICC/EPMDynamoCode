/** \file OuterCoreSource.hpp
 *  \brief Implementation of a field with a source term in outer core
 */

#ifndef OUTERCORESOURCE_HPP
#define OUTERCORESOURCE_HPP

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
    * \brief Implementation of a field with a source term in outer core
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSourceTraits Traits describing the field
    */
   template <typename TSimType, template <typename> class TSourceTraits> class OuterCoreSource: public TSourceTraits<TSimType>::DomainType
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
         OuterCoreSource(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~OuterCoreSource() {};

         /**
          * @brief Get values of outer core source term
          */
         const SpectralFieldType& ocSrc() const;

         /**
          * @brief Set values of outer core source term
          */
         SpectralFieldType&  rOcSrc();
         
      protected:

      private:
         /**
          * @brief Storage for Source term
          */
         SpectralFieldType mOcSrc;

   };

   template <typename TSimType, template <typename>  class TSourceTraits> OuterCoreSource<TSimType, TSourceTraits>::OuterCoreSource(SmartTruncation pTrunc, typename OuterCoreSource<TSimType,TSourceTraits>::TransformType &transform)
      : TSourceTraits<TSimType>::DomainType(pTrunc, transform), mOcSrc(pTrunc, true)
   {
   }

   template <typename TSimType, template <typename>  class TSourceTraits> const typename OuterCoreSource<TSimType, TSourceTraits>::SpectralFieldType& OuterCoreSource<TSimType, TSourceTraits>::ocSrc() const
   {
      return this->mOcSrc;
   }

   template <typename TSimType, template <typename>  class TSourceTraits> typename OuterCoreSource<TSimType, TSourceTraits>::SpectralFieldType& OuterCoreSource<TSimType, TSourceTraits>::rOcSrc()
   {
      return this->mOcSrc;
   }

}

#endif // OUTERCORESOURCE_HPP
