/** \file ETD1Iteration.hpp
 *  \brief Implemenation of the simple ETD1 iteration step. This is also the A step for ETD2RK method. 
 */

#ifndef ETD1ITERATION_HPP
#define ETD1ITERATION_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD1Iteration: public SchemeIteration<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for a smart pointer to ETDOperators
         typedef EPMSHARED_PTR<typename ETDSchemeTraits<TSimType>::Operators> SmartETDOperators;

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

   template <typename TSimType> ETD1Iteration<TSimType>::ETD1Iteration(EPMFloat nFactor, SmartETDOperators pOpM0, SmartETDOperators pOpM1)
      : mNFactor(nFactor), mpOpM0(pOpM0), mpOpM1(pOpM1)
   {
   }

   template <typename TSimType> void ETD1Iteration<TSimType>::compute(typename ETD1Iteration<TSimType>::ScalarType &rVar, typename ETD1Iteration<TSimType>::ScalarType &rNTerms)
   {
      // Get number of harmonic degrees
      int nL = rVar.trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = rVar.minL();

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {  
         rNTerms.rLShell(l) *= this->mNFactor;
         this->mpOpM1->multiplyOrders(rNTerms.rLShell(l), rNTerms.lshell(l),l);
         this->mpOpM0->affineOrders(rVar.rLShell(l), rVar.lshell(l), rNTerms.lshell(l), l);
         this->mpOpM0->extendOrders(rVar.rLShell(l), l);
      }
   }

}

#endif // ETD1ITERATION_HPP
