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
    * \tparam TBCType Type of the boundary condition
    */
   template <typename TOpType, BCType TBCType> class BoundedOperator: public BoundedOperatorBase<TOpType>
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
          */
         virtual void implementBCs(const Matrix& bcRows);

         /**
          * @brief Compute the operator the bounded operator based on input
          *
          * @param factor Multiplicative factor
          * @param polyOp Polynomial operator defining the operator
          *
          * \bug PolynomialOperator needs restructuring
          */
         virtual void constructBOperator(const DynamoFloat factor, const PolynomialOperator<POperator>& polyOp);

         /**
          * @brief Solve linear equation
          *
          * @param vector RHS of the linear equation
          */
         virtual void solve(Array& vector);

      protected:

      private:
         /**
          * @brief Extension matrix storage
          */
         ExtensionMatrix<TBCType>  mExt;
   };

   template <typename TOpType, BCType TBCType> BoundedOperator<TOpType, TBCType>::BoundedOperator(const int nBC, const int nTau, const int id)
      : BoundedOperatorBase<TOpType>(nTau, nBC, id), mExt()
   {
   }

   template <typename TOpType, BCType TBCType> inline void BoundedOperator<TOpType, TBCType>::implementBCs(const Matrix& bcRows)
   {
      // Build the extension matrices
      this->mExt.build(bcRows);
   }

   template <typename TOpType, BCType TBCType> void BoundedOperator<TOpType, TBCType>::constructBOperator(const DynamoFloat factor, const PolynomialOperator<POperator>& polyOp)
   {
      this->mExt.restrictOperator(this->rOp(), factor, polyOp.op());
   }

   template <typename TOpType, BCType TBCType> void BoundedOperator<TOpType, TBCType>::solve(Array& vector)
   {
      // Impose NH BC value if applicable
      this->mExt.prepareRHS(vector);

      // Call basic solve
      this->solveEquation(vector);

      // Extend solution to full truncation
      this->mExt.extend(vector);
   }

}

#endif // BOUNDEDOPERATOR_HPP
