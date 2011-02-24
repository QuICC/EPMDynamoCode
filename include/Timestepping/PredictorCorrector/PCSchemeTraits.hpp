/** \file PCSchemeTraits.hpp
 *  \brief Traits to describe the requested Predictor/Corrector scheme
 */

#ifndef PCSCHEMETRAITS_HPP
#define PCSCHEMETRAITS_HPP

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
   template <typename, template <typename> class> class IterativeScheme;
   template <typename> class ThetaMethodImplementation;
   template <typename> class ThetaInfluenceMethodImplementation;

   /**
    * \brief Traits to describe the request Predictor/Corrector scheme
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PCSchemeTraits
   {
      public:
         /// Typedef for the "general" setup traits
         typedef ThetaTraits<TSimType>  MethodTraits;

         /// Typedef for the simple timestepper (no influence matrix)
         typedef IterativeScheme<TSimType, ThetaMethodImplementation>   Timestepper;  

         /// Typedef for the timestepper with influence matrix step
         typedef IterativeScheme<TSimType, ThetaInfluenceMethodImplementation>   InfluenceTimestepper;  

         /// Typedef for the timestep control object type
         typedef PCTimestepControl<MethodTraits>   TimestepControl;
   };

}

#endif // PCSCHEMETRAITS_HPP
