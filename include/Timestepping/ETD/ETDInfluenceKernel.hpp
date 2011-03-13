/** \file ETDInfluenceKernel.hpp
 *  \brief Implemenation of influence of the kernel solution
 */

#ifndef ETDINFLUENCEKERNEL_HPP
#define ETDINFLUENCEKERNEL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/ETD/ETDSchemeTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of influence of the kernel solution
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETDInfluenceKernel
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
          * @param pOpM1 Pointer to operator M1
          * @param pOpM2 Pointer to operator M2
          */
         ETDInfluenceKernel(EPMFloat nFactor, SmartETDOperators pOpM1);

         /**
          * @brief Destructor
          */
         virtual ~ETDInfluenceKernel() {};

         /**
          * @brief Compute the influence of the kernel
          *
          * @param rKernel Kernel decomposition
          * @param l Harmonic degree
          */
         void computeInfluence(Array  &rKernel, const int l);
         
      protected:
         /**
          * @brief multiplicative factor required for the non linear terms
          */
         EPMFloat mNFactor;

         /**
          * @brief Pointer to the ETD operator M1
          */
         SmartETDOperators  mpOpM1;

      private:
   };

   template <typename TSimType> ETDInfluenceKernel<TSimType>::ETDInfluenceKernel(EPMFloat nFactor, SmartETDOperators pOpM1)
      : mNFactor(nFactor), mpOpM1(pOpM1)
   {
   }

   template <typename TSimType> void ETDInfluenceKernel<TSimType>::computeInfluence(Array &rKernel, const int l)
   {
      int rows = this->mpOpM1->harmOp(l).op().rows();
      rKernel.topRows(rows) = this->mNFactor*this->mpOpM1->harmOp(l).op() * rKernel.topRows(rows);
      this->mpOpM1->extendOrders(rKernel, l);
   }

}

#endif // ETDINFLUENCEKERNEL_HPP
