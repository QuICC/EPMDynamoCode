/** \file ETD2RKInfluenceKernel.hpp
 *  \brief Implemenation of influence of the kernel solution
 */

#ifndef ETD2RKINFLUENCEKERNEL_HPP
#define ETD2RKINFLUENCEKERNEL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/ETD/ETDIteration.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of influence of the kernel solution
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD2RKInfluenceKernel: public ETDIteration<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for a smart pointer to ETDOperators
         typedef EPMSHARED_PTR<typename ETDSchemeTraits<TSimType>::Operators> SmartETDOperators;

         /**
          * @brief Constructor
          *
          * @param pOpM1 Pointer to operator M1
          * @param pOpM2 Pointer to operator M2
          */
         ETD2RKInfluenceKernel(SmartETDOperators pOpM1, SmartETDOperators pOpM2);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKInfluenceKernel() {};

         /**
          * @brief Compute the influence of the kernel
          *
          * @param rKernel Kernel decomposition
          * @param l Harmonic degree
          */
         virtual void computeInfluence(Array  &rKernel, const int l);
         
      protected:
         /**
          * @brief Pointer to the ETD operator M1
          */
         SmartETDOperators  mpOpM1;

         /**
          * @brief Pointer to the ETD operator M2
          */
         SmartETDOperators  mpOpM2;

      private:
   };

   template <typename TSimType> ETD2RKInfluenceKernel<TSimType>::ETD2RKInfluenceKernel(SmartETDOperators pOpM1, SmartETDOperators pOpM2)
      : mpOpM1(pOpM1), mpOpM2(pOpM2)
   {
   }

   template <typename TSimType> void ETD2RKInfluenceKernel<TSimType>::computeInfluence(Array &rKernel, const int l)
   {
      rKernel = (this->mpOpM1->harmOp(l).op() + this->mpOpM2->harmOp(l).op()/this->h()) * rKernel;
   }

}

#endif // ETD2RKINFLUENCEKERNEL_HPP
