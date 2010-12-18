/** \file VelocityField.hpp
 *  \brief Implementation of the Magnetic vector field
 */

#ifndef VELOCITYFIELD_HPP
#define VELOCITYFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Traits/DefaultVelTraits.hpp"
#include "PhysicalFields/Traits/FieldNTermsTraits.hpp"
#include "PhysicalFields/PhysicalVariable.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Velocity vector field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TVelTraits Type of the field
    */
   template <typename TSimType, template <typename> class TVelTraits = DefaultVelTraits> class VelocityField: public PhysicalVariable<TSimType, TVelTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the type of the required non-linear terms
         typedef PhysicalVariable<TSimType,FieldNTermsTraits>    NTermsType;

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

   template <typename TSimType, template <typename> class TVelTraits> VelocityField<TSimType, TVelTraits>::VelocityField(SmartTruncation pTrunc, typename VelocityField<TSimType, TVelTraits>::TransformType &transform)
      : PhysicalVariable<TSimType, TVelTraits>(pTrunc, transform)
   {
   }

}

#endif // VELOCITYFIELD_HPP
