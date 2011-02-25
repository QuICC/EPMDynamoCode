/** \file ETD2RKTimestep.hpp
 *  \brief Implementation of the ETD2RK timestep step
 */

#ifndef ETD2RKTIMESTEP_HPP
#define ETD2RKTIMESTEP_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/SchemeIteration.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the ETD2RK timestep step
    *
    * This class performs the actual timestep comuputation after the intermediate
    * computations are performed.
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD2RKTimestep: public SchemeIteration<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// typedef for a pointer to a scalar type
         typedef EPMSHARED_PTR<ScalarType> SmartScalarType;

         /// Typedef for a smart pointer to ETDOperators
         typedef EPMSHARED_PTR<typename ETDSchemeTraits<TSimType>::Operators> SmartETDOperators;

         /**
          * @brief Constructor
          *
          * @param pOpM2 Pointer to the \f$M_2\f$ operator
          */
         ETD2RKTimestep(SmartScalarType pOldN, SmartETDOperators pOpM2);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKTimestep() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);
         
      protected:
         /**
          * @brief Storage for the old NTerms
          */
         SmartScalarType   mpOldN;

         /**
          * @brief ETD operator M2
          */
         SmartETDOperators  mpOpM2;

      private:
   };

   template <typename TSimType> ETD2RKTimestep<TSimType>::ETD2RKTimestep(typename ETD2RKTimestep<TSimType>::SmartScalarType pOldN, typename ETD2RKTimestep<TSimType>::SmartETDOperators pOpM2)
      : mpOldN(pOldN), mpOpM2(pOpM2)
   {
   }

   template <typename TSimType> void ETD2RKTimestep<TSimType>::compute(typename ETD2RKTimestep<TSimType>::ScalarType &rVar, typename ETD2RKTimestep<TSimType>::ScalarType &rNTerms)
   {
      // Get number of harmonic degrees
      int nL = rVar.trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = rVar.minL();

      EPMFloat h = this->mpOpM2->h();

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {
         rNTerms.rLShell(l) -= this->mpOldN->lshell(l);
         this->mpOpM2->multiplyOrders(rNTerms.rLShell(l), rNTerms.lshell(l), l);
         rVar.rLShell(l) += rNTerms.lshell(l)/h;
         this->mpOpM2->extendOrders(rNTerms.rLShell(l),l);
      }
   }

}

#endif // ETD2RKTIMESTEP_HPP
