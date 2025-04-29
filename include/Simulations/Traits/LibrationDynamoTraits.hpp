/** \file LibrationDynamoTraits.hpp
 *  \brief Traits for the precession dynamo simulation implementation
 */

#ifndef LIBRATIONDYNAMOTRAITS_HPP
#define LIBRATIONDYNAMOTRAITS_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

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

   class LibrationDynamoCodTraits
   {
      public:
         /// Typedef for the simulation type
         typedef SimulationConfig::NumericalScheme  SimulationType;

         /// Does the codensity field have an IC imposed field ? 
         static const bool HasICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasOCImposed = false; 
   };

   class LibrationDynamoMagTraits
   {
      public:
         /// Does the codensity field have an IC imposed field ? 
         static const bool HasICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasOCImposed = false; 
   };

   class LibrationDynamoVelTraits
   {
      public:
         /// Does the codensity field have an IC imposed field ? 
         static const bool HasICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasOCImposed = false; 
   };

   /**
    * @brief LibrationDynamo traits class
    */
   class LibrationDynamoTraits
   {
      public:
         /// Typedef for the codensity scalar type
         typedef CodensityScalar<>  CodType;

         /// Does the codensity field have an IC source ?
         static const bool HasCodICSource = false;

         /// Does the codensity field have an OC source ?
         static const bool HasCodOCSource = false;

         /// Does the codensity field have an IC imposed field ? 
         static const bool HasCodICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasCodOCImposed = false; 

         /// Does simulation required a codensity field
         static const bool NeedCodensity = false;

         /// Does simulation required a codensity gradient field
         static const bool NeedCodensityGrad = false;



         /// Typedef for the magnetic field type
         typedef MagneticField<>  MagType;

         /// Does the magnetic field have an IC source ?
         static const bool HasMagICSource = false;

         /// Does the magnetic field have an OC source ?
         static const bool HasMagOCSource = false;

         /// Does the magnetic field have an IC imposed field ? 
         static const bool HasMagICImposed = false; 

         /// Does the magnetic field have an OC imposed field ? 
         static const bool HasMagOCImposed = false; 

         /// Does simulation required a magnetic field
         static const bool NeedMagnetic = true;

         /// Does simulation required a magnetic curl field
         static const bool NeedMagneticCurl = true;



         /// Typedef for the velocity field type
         typedef VelocityField<>   VelType; 

         /// Does the velocity field have an IC source ?
         static const bool HasVelICSource = false;

         /// Does the velocity field have an OC source ?
         static const bool HasVelOCSource = false;

         /// Does the velocity field have an IC imposed field ? 
         static const bool HasVelICImposed = false; 

         /// Does the velocity field have an OC imposed field ? 
         static const bool HasVelOCImposed = false; 

         /// Does simulation required a velocity field
         static const bool NeedVelocity = true;

         /// Does simulation required a velocity curl field
         static const bool NeedVelocityCurl = true;
   };

}

#endif // LIBRATIONDYNAMOTRAITS_HPP
