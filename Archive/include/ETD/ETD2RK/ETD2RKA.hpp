/** \file ETD2RKA.hpp
 */

#ifndef ETD2RKA_HPP
#define ETD2RKA_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/ETD/ETDIteration.hpp"

namespace EPMDynamo {

   template <typename TSim> class ETD2RKA: public ETDIteration<TSim>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          */
         ETD2RKA(SmartETDOperator pOpM0, SmartETDOperator pOpM1);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKA();

         /**
          * @brief Compute the implemented iteration
          */
         virtual void compute(ScalarType &rVar, ScalarType rNTerms);
         
      protected:
         /**
          * @brief ETD operator M0
          */
         SmartETDOperator  mpOpM0;

         /**
          * @brief ETD operator M1
          */
         SmartETDOperator  mpOpM1;

      private:
   };

   template <typename TSim> ETD2RKA<TSim>::ETD2RKA(SmartETDOperator pOpM0, SmartETDOperator pOpM1)
      : mpOpM0(pOpM0), mpOpM1(pOpM1)
   {
   }

   template <typename TSim> ETD2RKA<TSim>::~ETD2RKA()
   {
   }

   template <typename TSim> void ETD2RKA<TSim>::compute(ETD2RK<TSim>::ScalarType &rVar, ETD2RK<TSim>::ScalarType &rNTerms)
   {
      // rVar = this->mpOpM0 * rVar + this->mpOpM1 * rNTerms;
   }

}

#endif // ETD2RKA_HPP
