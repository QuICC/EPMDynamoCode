/** \file DefaultCodDomainTraits.hpp
 *  \brief Traits for the Domain of the default codensity scalar field
 */

#ifndef DEFAULTCODDOMAINTRAITS_HPP
#define DEFAULTCODDOMAINTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/DefaultCodFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the Domain of the default codensity scalar field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class DefaultCodDomainTraits
   {
      public:
         /// Typedef for the spectral field type
         typedef typename DefaultCodFieldTraits<TSimType>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the domain type
         typedef FullSphereField<TSimType, DefaultCodFieldTraits> DomainType;
   };

}

#endif // DEFAULTCODDOMAINTRAITS_HPP
