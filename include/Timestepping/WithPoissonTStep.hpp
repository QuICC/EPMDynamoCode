/** \file WithPoissonTStep.hpp
 *  \brief Temporary implementation of a traits class...
 */

#ifndef WITHPOISSONTSTEP_HPP
#define WITHPOISSONTSTEP_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief Timestepping poisson traits
    *
    * \bug This has to be heavliy improved and restructured
    */
   template <typename T> class WithPoissonTStep
   {
      public:
         /// Define typedef to use the poisson step in timestepper
         typedef typename T::PoissonTimestepper   Timestepper;  
   };

}

#endif // WITHPOISSONTSTEP_HPP
