/** \file ThermalConvectionTraits.hpp
 *  \brief Traits for the dynamo simulation implementation
 */

#ifndef THERMALCONVECTIONTRAITS_HPP
#define THERMALCONVECTIONTRAITS_HPP

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
    * @brief Dynamo traits class
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThermalConvectionTraits
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

         /// Does simulation required a codensity gradient field
         static const bool NeedCodensityGrad = true;

         /// Does simulation required a magnetic field
         static const bool NeedMagnetic = false;

         /// Does simulation required a magnetic curl field
         static const bool NeedMagneticCurl = false;

         /// Does simulation required a velocity field
         static const bool NeedVelocity = true;

         /// Does simulation required a velocity curl field
         static const bool NeedVelocityCurl = true;
   };

}

#endif // THERMALCONVECTIONTRAITS_HPP
