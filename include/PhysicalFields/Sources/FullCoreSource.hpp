/** \file FullCoreSource.hpp
 *  \brief Implementation of a field with a source term in full core
 */

#ifndef FULLCORESOURCE_HPP
#define FULLCORESOURCE_HPP

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
    * \brief Implementation of a field with a source term in full core
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSourceTraits Traits describing the field
    */
   template <typename TSimType, template <typename> class TSourceTraits> class FullCoreSource: public TSourceTraits<TSimType>::DomainType
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

   template <typename TSimType, template <typename>  class TSourceTraits> FullCoreSource<TSimType, TSourceTraits>::FullCoreSource(SmartTruncation pTrunc, typename FullCoreSource<TSimType,TSourceTraits>::TransformType &transform)
      : TSourceTraits<TSimType>::DomainType(pTrunc, transform), mIcSrc(pTrunc, true), mOcSrc(pTrunc, true)
   {
   }

   template <typename TSimType, template <typename>  class TSourceTraits> const typename FullCoreSource<TSimType, TSourceTraits>::SpectralFieldType& FullCoreSource<TSimType, TSourceTraits>::icSrc() const
   {
      return this->mIcSrc;
   }

   template <typename TSimType, template <typename>  class TSourceTraits> typename FullCoreSource<TSimType, TSourceTraits>::SpectralFieldType& FullCoreSource<TSimType, TSourceTraits>::rIcSrc()
   {
      return this->mIcSrc;
   }

   template <typename TSimType, template <typename>  class TSourceTraits> const typename FullCoreSource<TSimType, TSourceTraits>::SpectralFieldType& FullCoreSource<TSimType, TSourceTraits>::ocSrc() const
   {
      return this->mOcSrc;
   }

   template <typename TSimType, template <typename>  class TSourceTraits> typename FullCoreSource<TSimType, TSourceTraits>::SpectralFieldType& FullCoreSource<TSimType, TSourceTraits>::rOcSrc()
   {
      return this->mOcSrc;
   }

   template <typename TSimType, template <typename>  class TSourceTraits> void FullCoreSource<TSimType, TSourceTraits>::initialiseZeros()
   {
      TSourceTraits<TSimType>::DomainType::initialiseZeros();

      this->mIcSrc.intialiseZeros();

      this->mOcSrc.intialiseZeros();
   }

}

#endif // FULLCORESOURCE_HPP
