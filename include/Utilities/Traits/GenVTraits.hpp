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
#include "PhysicalFields/Traits/GenCodTraits.hpp"
#include "PhysicalFields/Traits/GenMagTraits.hpp"
#include "PhysicalFields/Traits/GenVelTraits.hpp"

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
         typedef CodensityScalar<TSimType, GenCodTraits>  CodType;

         /// Does the codensity field have an IC source ?
         static const bool HasCodICSource = false;

         /// Does the codensity field have an OC source ?
         static const bool HasCodOCSource = true;

         /// Does the codensity field have an IC imposed field ? 
         static const bool HasCodICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasCodOCImposed = true; 

         /// Does simulation required a codensity field
         static const bool NeedCodensity = false;

         /// Does simulation required a codensity gradient field
         static const bool NeedCodensityGrad = false;



         /// Typedef for the magnetic field type
         typedef MagneticField<TSimType, GenMagTraits>  MagType;

         /// Does the magnetic field have an IC source ?
         static const bool HasMagICSource = false;

         /// Does the magnetic field have an OC source ?
         static const bool HasMagOCSource = true;

         /// Does the magnetic field have an IC imposed field ? 
         static const bool HasMagICImposed = false; 

         /// Does the magnetic field have an OC imposed field ? 
         static const bool HasMagOCImposed = true; 

         /// Does simulation required a magnetic field
         static const bool NeedMagnetic = false;

         /// Does simulation required a magnetic curl field
         static const bool NeedMagneticCurl = false;



         /// Typedef for the velocity field type
         typedef VelocityField<TSimType, GenVelTraits>   VelType; 

         /// Does the velocity field have an IC source ?
         static const bool HasVelICSource = false;

         /// Does the velocity field have an OC source ?
         static const bool HasVelOCSource = true;

         /// Does the velocity field have an IC imposed field ? 
         static const bool HasVelICImposed = false; 

         /// Does the velocity field have an OC imposed field ? 
         static const bool HasVelOCImposed = true; 

         /// Does simulation required a velocity field
         static const bool NeedVelocity = true;

         /// Does simulation required a velocity curl field
         static const bool NeedVelocityCurl = false;
   };

}

#endif // GENVTRAITS_HPP
