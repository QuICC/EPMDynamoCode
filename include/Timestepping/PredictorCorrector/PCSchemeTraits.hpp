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
   template <typename> class ThetaPoissonMethod;
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
         /// Typedef for the timestepper without poisson step
         typedef PCScheme<TSimType, ThetaMethod>   Timestepper;  

         /// Typedef for the timestepper with poisson step
         typedef PCScheme<TSimType, ThetaPoissonMethod>   PoissonTimestepper;  

         /// Typedef for the timestepper with poisson step
         typedef PCTimestepControl   TimestepControl;

         /// Typedef for the error norm type
         typedef PCL2MaxNorm<TSimType>   ErrorNormType;
   };

}

#endif // PCSCHEMETRAITS_HPP
