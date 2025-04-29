/** \file ETD1Iteration.hpp
 *  \brief Implemenation of the simple ETD1 iteration step. This is also the A step for ETD2RK method. 
 */

#ifndef ETD1ITERATION_HPP
#define ETD1ITERATION_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/SchemeIteration.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the simple ETD1 iteration step. This is also the A step for ETD2RK method. 
    */
   class ETD1Iteration: public SchemeIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for a smart pointer to ETDOperators
         typedef EPMSHARED_PTR<ETDSchemeTraits::Operators> SmartETDOperators;

         /**
          * @brief Constructor
          *
          * @param nFactor Non linear terms multiplicative factor
          * @param pOpM0 Pointer to operator M0
          * @param pOpM1 Pointer to operator M1
          */
         ETD1Iteration(EPMFloat nFactor, SmartETDOperators pOpM0, SmartETDOperators pOpM1);

         /**
          * @brief Destructor
          */
         virtual ~ETD1Iteration() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);
         
      protected:
         /**
          * @brief multiplicative factor required for the non linear terms
          */
         EPMFloat mNFactor;

         /**
          * @brief Pointer to the ETD operator M0
          */
         SmartETDOperators  mpOpM0;

         /**
          * @brief Pointer to the ETD operator M1
          */
         SmartETDOperators  mpOpM1;

      private:
   };

}

#endif // ETD1ITERATION_HPP
