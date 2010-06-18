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
#include "PhysicalFields/Extensions/SimpleField.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the abstract concept of a physical variable, independently of it being a vector field or a scalar field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    * \tparam TSource Special extension of the field
    */
   template <typename TSimType, template <typename> class TFieldType, template <typename, template <typename> class> class TSource = SimpleField> class PhysicalVariable : public TSource<TSimType, TFieldType>
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

         /**
          * @brief Get truncation information
          */
         const SmartTruncation trunc() const;
         
      protected:

      private:
   };

   template <typename TSimType, template <typename>  class TFieldType, template <typename, template <typename> class> class TSource> PhysicalVariable<TSimType, TFieldType, TSource>::PhysicalVariable(SmartTruncation pTrunc, typename PhysicalVariable<TSimType,TFieldType,TSource>::TransformType &transform)
      : TSource<TSimType, TFieldType>(pTrunc, transform)
   {
   }

   template <typename TSimType, template <typename>  class TFieldType, template <typename, template <typename> class> class TSource> const SmartTruncation PhysicalVariable<TSimType, TFieldType, TSource>::trunc() const
   {
      return this->oc().trunc();
   }

}

#endif // PHYSICALVARIABLE_HPP
