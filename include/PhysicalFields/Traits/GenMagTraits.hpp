/** \file GenMagTraits.hpp
 *  \brief Traits for the generators magnetic field
 */

#ifndef GENMAGTRAITS_HPP
#define GENMAGTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/GenMagDomainTraits.hpp"
#include "PhysicalFields/Sources/OuterCoreSource.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the generators magnetic field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class GenMagTraits
   {
      public:
         /// Typedef for the variable type
         typedef OuterCoreSource<TSimType, GenMagDomainTraits> VariableType;
   };

}

#endif // GENMAGTRAITS_HPP
