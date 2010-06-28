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
#include "Timestepping/ETD/ETDIteration.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the ETD2RK timestep step
    *
    * This class performs the actual timestep comuputation after the intermediate
    * computations are performed.
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD2RKTimestep: public ETDIteration<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param pNn Pointer to the scalar
          * @param pOpM2 Pointer to the \f$M_2\f$ operator
          */
         ETD2RKTimestep(SmartScalarType pNn, SmartETDOperator pOpM2);

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

   template <typename TSimType> ETD2RKTimestep<TSimType>::ETD2RKTimestep(SmartScalarType pNn, SmartETDOperator pOpM2)
      : mpNTermsN(pNn), mpOpM2(pOpM2)
   {
   }

   template <typename TSimType> void ETD2RKTimestep<TSimType>::compute(ETD2RK<TSimType>::ScalarType &rVar, ETD2RK<TSimType>::ScalarType &rNTerms)
   {
      // rVar = rVar + this->mpOpM2 * (rNTerms - this->mpNTermsN)/h;
   }

}

#endif // ETD2RKTIMESTEP_HPP
