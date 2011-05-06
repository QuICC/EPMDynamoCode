/** \file PCSchemeTraits.hpp
 *  \brief Traits to describe the requested Predictor/Corrector scheme
 */

#ifndef PCSCHEMETRAITS_HPP
#define PCSCHEMETRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/PredictorCorrector/PCTimestepControl.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaTraits.hpp"

namespace EPMDynamo {

   // Forward declarations
   template <typename> class IterativeScheme;

   /**
    * \brief Traits to describe the request Predictor/Corrector scheme
    */
   class PCSchemeTraits
   {
      public:
         /// Typedef for the "general" setup traits
         typedef ThetaTraits  MethodTraits;

         /// Typedef for the simple timestepper (no influence matrix)
         typedef IterativeScheme<ThetaMethod>   Timestepper;  

         /// Typedef for the timestepper with influence matrix step
         typedef IterativeScheme<ThetaInfluenceMethod>   InfluenceTimestepper;  

         /// Typedef for the timestep control object type
         typedef PCTimestepControl<MethodTraits>   TimestepControl;
   };

}

#endif // PCSCHEMETRAITS_HPP
