/** \file ETDSchemeTraits.hpp
 *  \brief Traits to describe the requested ETD scheme
 */

#ifndef ETDSCHEMETRAITS_HPP
#define ETDSCHEMETRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Operators/DenseOperator.hpp"
#include "Timestepping/ETD/ETDTimestepControl.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKTraits.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"

namespace EPMDynamo {

   // Forward declaration
   template <typename, template <typename> class> class ETDScheme;
   template <typename> class ETD2RKMethod;
   template <typename> class ETD2RKInfluenceMethod;

   /**
    * \brief Traits to describe the requested ETD scheme
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETDSchemeTraits
   {
      public:
         /// Typedef for the "general" setup traits
         typedef ETD2RKTraits<TSimType>  MethodTraits;

         /// Typedef for the simple timestepper (without influence matrix step)
         typedef ETDScheme<TSimType, ETD2RKMethod>   Timestepper;  

         /// Typedef for the timestepper with influence matrix step
         typedef ETDScheme<TSimType, ETD2RKInfluenceMethod>   InfluenceTimestepper;  

         /// Typedef for the timestep control object
         typedef ETDTimestepControl<MethodTraits>   TimestepControl;

         /// Typedefs for the ETD operators
         typedef ETDOperators<TSimType, DenseOperator>  Operators;
   };

}

#endif // ETDSCHEMETRAITS_HPP
