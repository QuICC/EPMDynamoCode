/** \file ETDSchemeTraits.hpp
 *  \brief Traits to describe the requested ETD scheme
 */

#ifndef ETDSCHEMETRAITS_HPP
#define ETDSCHEMETRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Operators/DenseOperator.hpp"
#include "Timestepping/ETD/ETDTimestepControl.hpp"
#include "Timestepping/ETD/ETD1/ETD1Traits.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKTraits.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits to describe the requested ETD scheme
    */
   class ETDSchemeTraits
   {
      public:
         /// Typedef for the "general" setup traits
//         typedef ETD1Traits  MethodTraits;
         typedef ETD2RKTraits  MethodTraits;

         /// Typedef for the simple timestepper (without influence matrix step)
//         typedef IterativeScheme<ETD1Method>   Timestepper;  
         typedef IterativeScheme<ETD2RKMethod>   Timestepper;  

         /// Typedef for the timestepper with influence matrix step
//         typedef IterativeScheme<ETD1InfluenceMethod>   InfluenceTimestepper;  
         typedef IterativeScheme<ETD2RKInfluenceMethod>   InfluenceTimestepper;  

         /// Typedef for the timestep control object
         typedef ETDTimestepControl<MethodTraits>   TimestepControl;

         /// Typedefs for the ETD operators
         typedef ETDOperators<DenseOperator>  Operators;
   };

}

#endif // ETDSCHEMETRAITS_HPP
