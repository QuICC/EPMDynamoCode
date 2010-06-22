/** \file GenCTraits.hpp
 *  \brief Traits for a codensity generator
 */

#ifndef GENCTRAITS_HPP
#define GENCTRAITS_HPP

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
    * \brief Traits for a codensity generator
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class GenCTraits
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
         static const bool NeedMagnetic = false;

         /// Does simulation required a velocity field
         static const bool NeedVelocity = false;
   };

}

#endif // GENCTRAITS_HPP
