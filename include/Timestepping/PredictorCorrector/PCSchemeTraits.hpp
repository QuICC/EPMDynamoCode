/** \file PCSchemeTraits.hpp
 *  \brief Traits to describe the request Predictor/Corrector scheme
 */

#ifndef PCSCHEMETRAITS_HPP
#define PCSCHEMETRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   // Forward declarations
   class PCTimestepControl;
   template <typename, template <typename> class> class PCScheme;
   template <typename> class ThetaMethod;
   template <typename> class ThetaInfluenceMethod;
   template <typename> class ErrorL2Max;

   /**
    * \brief Traits to describe the request Predictor/Corrector scheme
    *
    * \tparam TSimType Type of the simulation
    *
    * \bug Needs a big cleaning and restructuration
    */
   template <typename TSimType> class PCSchemeTraits
   {
      public:
         /// Typedef for the simple timestepper (no influence matrix)
         typedef PCScheme<TSimType, ThetaMethod>   Timestepper;  

         /// Typedef for the timestepper with influence matrix step
         typedef PCScheme<TSimType, ThetaInfluenceMethod>   InfluenceTimestepper;  

         /// Typedef for the timestep control object type
         typedef PCTimestepControl   TimestepControl;

         /// Typedef for the error norm type
         typedef ErrorL2Max<TSimType>   ErrorNormType;

         /// Type of timestep controller to use
         static const TimestepCtrlTypes  CtrlType = ElementaryCtrl;

         /// Order of the timestep scheme
         static const int order = Timestepper::order;
   };

}

#endif // PCSCHEMETRAITS_HPP
