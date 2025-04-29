/** \file DenseLUOperator.hpp
 *  \brief Implementation of a dense LU factorised (with full pivoting) operator
 */

#ifndef DENSELUOPERATOR_HPP
#define DENSELUOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Operators/DenseOperator.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a dense LU factorised (with full pivoting) operator
    */
   class DenseLUOperator: public DenseOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size Size of the square operator
          */
         DenseLUOperator(const int size);

         /**
          * @brief Destructor
          */
         virtual ~DenseLUOperator() {};

         /**
          * @brief Finalise operator computation
          */
         virtual void finaliseOp();
         
      protected:
         /**
          * @brief Get pivoting indexes
          */
         const ArrayI&  ipiv() const;

         /**
          * @brief Set the pivoting indexes
          */
         ArrayI&  rIpiv();

         /**
          * @brief Compute LU factorisation
          */
         void factorise();

         /**
          * @brief Solve linear equation
          *
          * @param vector RHS of linear equation
          */
         void solveEquation(Array& vector);

         /**
          * @brief Compute inverse from LU factorisation
          */
         void inverse();

      private:
         /**
          * @brief Pivoting indexes storage
          */
         ArrayI   mIpiv;
   };

   inline const ArrayI& DenseLUOperator::ipiv() const
   {
      return this->mIpiv;
   }

   inline ArrayI& DenseLUOperator::rIpiv()
   {
      return this->mIpiv;
   }

   inline void DenseLUOperator::finaliseOp()
   {
      this->factorise();
   }

}

#endif // DENSELUOPERATOR_HPP
