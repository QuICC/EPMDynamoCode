/** \file WithoutPoissonTStep.hpp
 *  \brief Temporary implementation of a traits class...
 */

#ifndef WITHOUTPOISSONTSTEP_HPP
#define WITHOUTPOISSONTSTEP_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief Timestepping without poisson step traits
    *
    * \bug This has to be heavliy improved and restructured
    */
   template <typename T> class WithoutPoissonTStep
   {
      public:
         /// Define typedef for standard timestepper without poisson
         typedef typename T::Timestepper   Timestepper;  
   };

}

#endif // WITHOUTPOISSONTSTEP_HPP
