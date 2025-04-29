/** \file VelocityField.hpp
 *  \brief Implementation of the Magnetic vector field
 */

#ifndef VELOCITYFIELD_HPP
#define VELOCITYFIELD_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/DefaultVelTraits.hpp"
#include "PhysicalFields/Traits/FieldNTermsTraits.hpp"
#include "PhysicalFields/PhysicalVariable.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Velocity vector field
    *
    * \tparam TVelTraits Type of the field
    */
   template <typename TVelTraits = DefaultVelTraits> class VelocityField: public PhysicalVariable<TVelTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the type of the required non-linear terms
         typedef PhysicalVariable<FieldNTermsTraits>    NTermsType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         VelocityField(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VelocityField() {};
         
      protected:

      private:
   };

   template <typename TVelTraits> VelocityField<TVelTraits>::VelocityField(SmartTruncation pTrunc, typename VelocityField<TVelTraits>::TransformType &transform)
      : PhysicalVariable<TVelTraits>(pTrunc, transform)
   {
   }

}

#endif // VELOCITYFIELD_HPP
