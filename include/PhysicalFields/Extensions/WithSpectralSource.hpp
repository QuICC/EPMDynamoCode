/** \file WithSpectralSource.hpp
 *  \brief Implementation of a field with spectral source
 */

#ifndef WITHSPECTRALSOURCE_HPP
#define WITHSPECTRALSOURCE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"
#include "PhysicalFields/Extensions/SimpleField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a field with spectral source
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    */
   template <typename TSimType, template <typename> class TFieldType> class WithSpectralSource: public SimpleField<TSimType, TFieldType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Transform object
          */
         WithSpectralSource(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~WithSpectralSource() {};

         /**
          * @brief Get values of source term
          */
         const typename TFieldType<TSimType>::SpectralFieldType& src() const;

         /**
          * @brief Set values of source term
          */
         typename TFieldType<TSimType>::SpectralFieldType&  rSrc();
         
      protected:

      private:
         /**
          * @brief Storage for Source term
          */
         typename TFieldType<TSimType>::SpectralFieldType mSrc;

   };

   template <typename TSimType, template <typename>  class TFieldType> WithSpectralSource<TSimType, TFieldType>::WithSpectralSource(SmartTruncation pTrunc, typename WithSpectralSource<TSimType,TFieldType>::TransformType &transform)
      : SimpleField<TSimType, TFieldType>(pTrunc, transform), mSrc(pTrunc)
   {
   }

   template <typename TSimType, template <typename>  class TFieldType> const typename TFieldType<TSimType>::SpectralFieldType& WithSpectralSource<TSimType, TFieldType>::src() const
   {
      return this->mSrc;
   }

   template <typename TSimType, template <typename>  class TFieldType> typename TFieldType<TSimType>::SpectralFieldType& WithSpectralSource<TSimType, TFieldType>::rSrc()
   {
      return this->mSrc;
   }

}

#endif // WITHSPECTRALSOURCE_HPP
