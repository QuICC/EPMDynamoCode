/** \file PhysicalVariable.hpp
 *  \brief Implementation of a general physical variable
 */

#ifndef PHYSICALVARIABLE_HPP
#define PHYSICALVARIABLE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the abstract concept of a physical variable, independently of it being a vector field or a scalar field
    *
    * \tparam TVariableTraits Traits describing properties of field (domain, sources, imposed, ...)
    */
   template <typename TVariableTraits> class PhysicalVariable : public TVariableTraits::VariableType
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         PhysicalVariable(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Simple empty destructor
          */
         virtual ~PhysicalVariable() {};
         
      protected:

      private:
   };

   template <typename TVariableTraits> PhysicalVariable<TVariableTraits>::PhysicalVariable(SmartTruncation pTrunc, typename PhysicalVariable<TVariableTraits>::TransformType &transform)
      : TVariableTraits::VariableType(pTrunc, transform)
   {
   }

}

#endif // PHYSICALVARIABLE_HPP
