/** \file OuterCoreSource.hpp
 *  \brief Implementation of a field with a source term in outer core
 */

#ifndef OUTERCORESOURCE_HPP
#define OUTERCORESOURCE_HPP

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
    * \brief Implementation of a field with a source term in outer core
    *
    * \tparam TSourceTraits Traits describing the field
    */
   template <typename TSourceTraits> class OuterCoreSource: public TSourceTraits::DomainType
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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:

      private:
         /**
          * @brief Storage for Source term
          */
         SpectralFieldType mOcSrc;

   };

   template <typename TSourceTraits> OuterCoreSource<TSourceTraits>::OuterCoreSource(SmartTruncation pTrunc, typename OuterCoreSource<TSourceTraits>::TransformType &transform)
      : TSourceTraits::DomainType(pTrunc, transform), mOcSrc(pTrunc, true)
   {
   }

   template <typename TSourceTraits> const typename OuterCoreSource<TSourceTraits>::SpectralFieldType& OuterCoreSource<TSourceTraits>::ocSrc() const
   {
      return this->mOcSrc;
   }

   template <typename TSourceTraits> typename OuterCoreSource<TSourceTraits>::SpectralFieldType& OuterCoreSource<TSourceTraits>::rOcSrc()
   {
      return this->mOcSrc;
   }

   template <typename TSourceTraits> void OuterCoreSource<TSourceTraits>::initialiseZeros()
   {
      TSourceTraits::DomainType::initialiseZeros();

      this->mOcSrc.initialiseZeros();
   }

}

#endif // OUTERCORESOURCE_HPP
