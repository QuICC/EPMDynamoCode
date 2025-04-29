/** \file DenseOperator.hpp
 *  \brief General implementation of a dense operator
 */

#ifndef DENSEOPERATOR_HPP
#define DENSEOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a dense operator
    */
   class DenseOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows Number of rows of operator
          * @param cols Number of columns of operator
          */
         DenseOperator(const int rows, const int cols);

         /**
          * @brief Constructor
          *
          * @param size Size of a square operator
          * @param isFactorised Flag to trigger factorisation of operator
          */
         DenseOperator(const int size, const bool isFactorised = false);

         /**
          * @brief Destructor
          */
         virtual ~DenseOperator() {};

         /**
          * @brief Get operator
          */
         const Matrix&  op() const;
         
         /**
          * @brief Set operator
          */
         Matrix&  rOp();

         /**
          * @brief Finalise operator computation
          *
          * This method does nothing but is required for a more general implementation
          */
         virtual void  finaliseOp() {};

      protected:
         /**
          * @brief Get factorisation flag
          */
         bool isFactorised() const;

         /**
          * @brief Solve linear equation
          *
          * @param vector RHS of linear equation
          */
         void solveEquation(Array& vector);

         /**
          * @brief Invert the operator matrix
          */
         void inverse();

      private:
         /**
          * @brief flag for factorisation
          */
         bool mIsFactorised;

         /**
          * @brief Operator storage
          */
         Matrix   mOperator;
   };

   inline bool DenseOperator::isFactorised() const
   {
      return this->mIsFactorised;
   }

   inline const Matrix& DenseOperator::op() const
   {
      return this->mOperator;
   }

   inline Matrix& DenseOperator::rOp()
   {
      return this->mOperator;
   }

}

#endif // DENSEOPERATOR_HPP
