/** \file GenCodDomainTraits.hpp
 *  \brief Traits for the Domain of the generators codensity scalar field
 */

#ifndef GENCODDOMAINTRAITS_HPP
#define GENCODDOMAINTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/GenCodFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the Domain of the generators codensity scalar field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class GenCodDomainTraits
   {
      public:
         /// Typedef for the spectral field type
         typedef typename GenCodFieldTraits<TSimType>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the domain type
         typedef FullSphereField<TSimType, GenCodFieldTraits> DomainType;
   };

}

#endif // GENCODDOMAINTRAITS_HPP
