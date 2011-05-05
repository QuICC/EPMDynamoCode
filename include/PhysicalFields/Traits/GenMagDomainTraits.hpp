/** \file GenMagDomainTraits.hpp
 *  \brief Traits for the generators magnetic field domain
 */

#ifndef GENMAGDOMAINTRAITS_HPP
#define GENMAGDOMAINTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/GenMagFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the generators magnetic field domain
    */
   class GenMagDomainTraits
   {
      public:
         /// Typedef for the spectral field type
         typedef typename GenMagFieldTraits::SpectralFieldType  SpectralFieldType;

         /// Typedef for the variable type
         typedef FullSphereField<GenMagFieldTraits> DomainType;
   };

}

#endif // GENMAGDOMAINTRAITS_HPP
