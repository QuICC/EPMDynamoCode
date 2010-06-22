/** \file GenVTraits.hpp
 *  \brief Traits for a velocity generator
 */

#ifndef GENVTRAITS_HPP
#define GENVTRAITS_HPP

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
    * \brief Traits for a velocity generator
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class GenVTraits
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
         static const bool NeedMagnetic = false;

         /// Does simulation required a velocity field
         static const bool NeedVelocity = true;
   };

}

#endif // GENVTRAITS_HPP
