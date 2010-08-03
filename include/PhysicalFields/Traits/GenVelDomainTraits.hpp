/** \file GenVelDomainTraits.hpp
 *  \brief Traits for the generators velocity field domain
 */

#ifndef GENVELDOMAINTRAITS_HPP
#define GENVELDOMAINTRAITS_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class GenVelDomainTraits
   {
      public:
         /// Typedef for the spectral field type
         typedef typename GenVelFieldTraits<TSimType>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the variable type
         typedef FullSphereField<TSimType, GenVelFieldTraits> DomainType;
   };

}

#endif // GENVELDOMAINTRAITS_HPP
