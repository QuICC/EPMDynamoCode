/** \file GenVelDomainTraits.hpp
 *  \brief Traits for the generators velocity field domain
 */

#ifndef GENVELDOMAINTRAITS_HPP
#define GENVELDOMAINTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/GenVelFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the generators velocity field domain
    */
   class GenVelDomainTraits
   {
      public:
         /// Typedef for the spectral field type
         typedef typename GenVelFieldTraits::SpectralFieldType  SpectralFieldType;

         /// Typedef for the variable type
         typedef FullSphereField<GenVelFieldTraits> DomainType;
   };

}

#endif // GENVELDOMAINTRAITS_HPP
