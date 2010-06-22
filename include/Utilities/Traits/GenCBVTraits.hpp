/** \file GenCBVTraits.hpp
 *  \brief Traits for codensity, magnetic and velocity generator
 */

#ifndef GENCBVTRAITS_HPP
#define GENCBVTRAITS_HPP

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
    * \brief Traits for codensity, magnetic and velocity generator
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class GenCBVTraits
   {
      public:
         /// Typedef for the codensity scalar type
         typedef CodensityScalar<TSimType>  CodType;

         /// Typedef for the magnetic field type
         typedef MagneticField<TSimType>  MagType;

         /// Typedef for the velocity field type
         typedef VelocityField<TSimType>   VelType; 

         /// Does simulation required a codensity field
         static const bool NeedCodensity = true;

         /// Does simulation required a magnetic field
         static const bool NeedMagnetic = true;

         /// Does simulation required a velocity field
         static const bool NeedVelocity = true;
   };

}

#endif // GENCBVTRAITS_HPP
