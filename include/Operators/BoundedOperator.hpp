/** \file BoundedOperator.hpp
 *  \brief Implementation of a general bounded operator
 */

#ifndef BOUNDEDOPERATOR_HPP
#define BOUNDEDOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Operators/BoundedOperatorBase.hpp"
#include "Operators/ExtensionMatrix.hpp"
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/PolyOps/POperator.hpp"
#include "LAPACK_Iface.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a general bounded operator (homogeneous or non homogeneous)
    *
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class BoundedOperator: public BoundedOperatorBase<TOpType>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nBC Number of boundary conditions
          * @param nTau Size of the bounded operator
          * @param id ID of the operator (for example harmonic degree l)
          */
         BoundedOperator(const int nBC, const int nTau, const int id);

         /**
          * @brief Destructor
          */
         virtual ~BoundedOperator() {};

         /**
          * @brief Implement the boundary conditions
          *
          * @param bcRows Matrix of the boundary imposing row values
          * @param needNHMatrix Needs the NH imposing matrix ?
          */
         virtual void implementBCs(const Matrix& bcRows, bool needNHMatrix);

         /**
          * @brief Compute the operator the bounded operator based on input
          *
          * @param factor Multiplicative factor
          * @param polyOp Polynomial operator defining the operator
          */
         virtual void constructBOperator(const EPMFloat factor, const PolynomialOperator<POperator>& polyOp);

         /**
          * @brief Solve linear equation
          *
          * @param vector RHS of the linear equation
          * @param isHomogeneous Is the boundary condition homogeneous?
          */
         virtual void solve(Array& vector, bool isHomogeneous);

         /**
          * @brief Solve linear equation with zero BC value
          *
          * @param vector RHS of the linear equation
          */
         virtual void solveZero(Array& vector);

      protected:

      private:
         /**
          * @brief Extension matrix storage
          */
         ExtensionMatrix  mExt;

         /**
          * @brief Flag to test if NH matrix is needed
          */
         bool mNeedNHMatrix;
   };

   template <typename TOpType> BoundedOperator<TOpType>::BoundedOperator(const int nBC, const int nTau, const int id)
      : BoundedOperatorBase<TOpType>(nBC, nTau, id), mExt()
   {
   }

   template <typename TOpType> inline void BoundedOperator<TOpType>::implementBCs(const Matrix& bcRows, bool needNHMatrix)
   {
      // Build the extension matrices
      this->mExt.build(bcRows);

      // Set the NH matrix construction flag
      this->mNeedNHMatrix = needNHMatrix;
   }

   template <typename TOpType> void BoundedOperator<TOpType>::constructBOperator(const EPMFloat factor, const PolynomialOperator<POperator>& polyOp)
   {
      this->mExt.restrictOperator(this->rOp(), factor, polyOp.op());

      if(this->mNeedNHMatrix)
      {
         this->mExt.buildNHMatrix(factor, polyOp.op());
      }
   }

   template <typename TOpType> void BoundedOperator<TOpType>::solveZero(Array& vector)
   {
      // Call basic solve
      this->solveEquation(vector);

      // Extend solution to full truncation
      this->mExt.extendZero(vector);
   }

   template <typename TOpType> void BoundedOperator<TOpType>::solve(Array& vector, bool isHomogeneous)
   {
      if(isHomogeneous)
      {
         // Call basic solve
         this->solveEquation(vector);

         // Extend solution to full truncation
         this->mExt.extendZero(vector);
      } else
      {
         // Impose NH BC value if applicable
         this->mExt.prepareRHS(vector);

         // Call basic solve
         this->solveEquation(vector);

         // Extend solution to full truncation
         this->mExt.extend(vector);
      }
   }

}

#endif // BOUNDEDOPERATOR_HPP
