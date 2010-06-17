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
#include "Simulations/SimulationTraits.hpp"
#include "PhysicalFields/PhysicalNoDivField.hpp"
#include "PhysicalFields/PhysicalVariable.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Velocity vector field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    */
   template <typename TSimType, template <typename> class TFieldType = PhysicalNoDivField> class VelocityField: public PhysicalVariable<TSimType, TFieldType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

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

   template <typename TSimType, template <typename> class TFieldType> VelocityField<TSimType, TFieldType>::VelocityField(SmartTruncation pTrunc, typename VelocityField<TSimType, TFieldType>::TransformType &transform)
      : PhysicalVariable<TSimType, TFieldType>(pTrunc, transform)
   {
   }

}

#endif // VELOCITYFIELD_HPP
