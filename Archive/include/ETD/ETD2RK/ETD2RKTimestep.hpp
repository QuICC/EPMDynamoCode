/** \file ETD2RKTimestep.hpp
 */

#ifndef ETD2RKTIMESTEP_HPP
#define ETD2RKTIMESTEP_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/ETD/ETDIteration.hpp"

namespace EPMDynamo {

   template <typename TSim> class ETD2RKTimestep: public ETDIteration<TSim>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          */
         ETD2RKTimestep(SmartScalarType pNn, SmartETDOperator pOpM2);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKTimestep();

         /**
          * @brief Compute the implemented iteration
          */
         virtual void compute(ScalarType &rVar, ScalarType rNTerms);
         
      protected:
         /**
          * @brief Storage for the old NTerms
          */
         SmartScalarType   mpNTermsN;

         /**
          * @brief ETD operator M2
          */
         SmartETDOperator  mpOpM2;

      private:
   };

   template <typename TSim> ETD2RKTimestep<TSim>::ETD2RKTimestep(SmartScalarType pNn, SmartETDOperator pOpM2)
      : mpNTermsN(pNn), mpOpM2(pOpM2)
   {
   }

   template <typename TSim> ETD2RKTimestep<TSim>::~ETD2RKTimestep()
   {
   }

   template <typename TSim> void ETD2RKTimestep<TSim>::compute(ETD2RK<TSim>::ScalarType &rVar, ETD2RK<TSim>::ScalarType &rNTerms)
   {
      // rVar = rVar + this->mpOpM2 * (rNTerms - this->mpNTermsN)/h;
   }

}

#endif // ETD2RKTIMESTEP_HPP
