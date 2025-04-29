/** \file DynamoTraits.hpp
 *  \brief Traits for the dynamo simulation implementation
 */

#ifndef DYNAMOTRAITS_HPP
#define DYNAMOTRAITS_HPP

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

   class DynamoCodTraits
   {
      public:
         /// Typedef for the simulation type
         typedef SimulationConfig::NumericalScheme  SimulationType;

         /// Does the codensity field have an IC imposed field ? 
         static const bool HasICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasOCImposed = false; 
   };

   class DynamoMagTraits
   {
      public:
         /// Does the codensity field have an IC imposed field ? 
         static const bool HasICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasOCImposed = false; 
   };

   class DynamoVelTraits
   {
      public:
         /// Does the codensity field have an IC imposed field ? 
         static const bool HasICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasOCImposed = false; 
   };

   /**
    * @brief Dynamo traits class
    */
   class DynamoTraits
   {
      public:
         /// Typedef for the codensity scalar type
         typedef CodensityScalar<>  CodType;

         /// Does the codensity field have an IC source ?
         static const bool HasCodICSource = false;

         /// Does the codensity field have an OC source ?
         static const bool HasCodOCSource = true;

         /// Does the codensity field have an IC imposed field ? 
         static const bool HasCodICImposed = false; 

         /// Does the codensity field have an OC imposed field ? 
         static const bool HasCodOCImposed = false; 

         /// Does simulation required a codensity field
         static const bool NeedCodensity = true;

         /// Does simulation required a codensity gradient field
         static const bool NeedCodensityGrad = true;



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

#endif // DYNAMOTRAITS_HPP
