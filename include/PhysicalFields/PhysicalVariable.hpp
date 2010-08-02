/** \file PhysicalVariable.hpp
 *  \brief Implementation of a general physical variable
 */

#ifndef PHYSICALVARIABLE_HPP
#define PHYSICALVARIABLE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"

#include "boost/static_assert.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the abstract concept of a physical variable, independently of it being a vector field or a scalar field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TVariableTraits Traits describing properties of field (domain, sources, imposed, ...)
    */
   template <typename TSimType, template <typename> class TVariableTraits> class PhysicalVariable : public TVariableTraits<TSimType>::VariableType
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
         PhysicalVariable(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Simple empty destructor
          */
         virtual ~PhysicalVariable() {};
         
      protected:

      private:
   };

   template <typename TSimType, template <typename>  class TVariableTraits> PhysicalVariable<TSimType, TVariableTraits>::PhysicalVariable(SmartTruncation pTrunc, typename PhysicalVariable<TSimType,TVariableTraits>::TransformType &transform)
      : TVariableTraits<TSimType>::VariableType(pTrunc, transform)
   {
   }

}

#endif // PHYSICALVARIABLE_HPP
