/** \file GenBTraits.hpp
 *  \brief Traits for a magnetic generator
 */

#ifndef GENBTRAITS_HPP
#define GENBTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "PhysicalFields/CodensityScalar.hpp"
#include "PhysicalFields/MagneticField.hpp"
#include "PhysicalFields/VelocityField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for a magnetic generator
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class GenBTraits
   {
      public:
         /// Typedef for the codensity scalar type
         typedef CodensityScalar<TSimType>  CodType;

         /// Typedef for the magnetic field type
         typedef MagneticField<TSimType>  MagType;

         /// Typedef for the velocity field type
         typedef VelocityField<TSimType>   VelType; 

         /// Does simulation required a codensity field
         static const bool NeedCodensity = false;

         /// Does simulation required a magnetic field
         static const bool NeedMagnetic = true;

         /// Does simulation required a velocity field
         static const bool NeedVelocity = false;
   };

}

#endif // GENBTRAITS_HPP
