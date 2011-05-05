/** \file DefaultCodDomainTraits.hpp
 *  \brief Traits for the Domain of the default codensity scalar field
 */

#ifndef DEFAULTCODDOMAINTRAITS_HPP
#define DEFAULTCODDOMAINTRAITS_HPP

// Configuration includes
//

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
    */
   class DefaultCodDomainTraits
   {
      public:
         /// Typedef for the spectral field type
         typedef DefaultCodFieldTraits::SpectralFieldType  SpectralFieldType;

         /// Typedef for the domain type
         typedef FullSphereField<DefaultCodFieldTraits> DomainType;
   };

}

#endif // DEFAULTCODDOMAINTRAITS_HPP
