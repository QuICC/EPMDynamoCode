/** \file CodensityDiffusionTraits.hpp
 *  \brief Traits for the codensity diffusion simulation implementation
 */

#ifndef CODENSITYDIFFUSIONTRAITS_HPP
#define CODENSITYDIFFUSIONTRAITS_HPP

// Configuration includes
//

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
    * \brief Traits for the codensity diffusion simulation implementation
    */
   class CodensityDiffusionTraits
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
         static const bool NeedMagnetic = false;

         /// Does simulation required a magnetic curl field
         static const bool NeedMagneticCurl = false;



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
         static const bool NeedVelocity = false;

         /// Does simulation required a velocity curl field
         static const bool NeedVelocityCurl = false;
   };

}

#endif // CODENSITYDIFFUSIONTRAITS_HPP
