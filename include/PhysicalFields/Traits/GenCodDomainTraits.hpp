/** \file GenCodDomainTraits.hpp
 *  \brief Traits for the Domain of the generators codensity scalar field
 */

#ifndef GENCODDOMAINTRAITS_HPP
#define GENCODDOMAINTRAITS_HPP

// Configuration includes
//

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
    */
   class GenCodDomainTraits
   {
      public:
         /// Typedef for the spectral field type
         typedef GenCodFieldTraits::SpectralFieldType  SpectralFieldType;

         /// Typedef for the domain type
         typedef FullSphereField<GenCodFieldTraits> DomainType;
   };

}

#endif // GENCODDOMAINTRAITS_HPP
