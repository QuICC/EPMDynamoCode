/** \file BCRowOperator.hpp
 *  \brief Implementation of a boundary bordering bounded operator
 */

#ifndef BCROWOPERATOR_HPP
#define BCROWOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Operators/BoundedOperatorBase.hpp"
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/PolyOps/POperator.hpp"
#include "LAPACK_Iface.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a boundary bordering bounded operator
    *
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class BCRowOperator: public BoundedOperatorBase<TOpType>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nBC Number of boundary conditions
          * @param nTau Size of the bounded operator
          * @param id ID of the operator (for example harmonic degree l)
          */
         BCRowOperator(const int nBC, const int nTau, const int id);

         /**
          * @brief Destructor
          */
         virtual ~BCRowOperator() {};

         /**
          * @brief Implement the boundary conditions
          *
          * @param bcRows Matrix of the boundary imposing row values
          */
         virtual void implementBCs(const Matrix& bcRows);

         /**
          * @brief Set operator's diagonal value
          *
          * @param factor Multiplicative factor
          */
         virtual void setOperator(const EPMFloat factor);

         /**
          * @brief Set operator value
          *
          * @param factor Multiplicative factor
          * @param mat    Input operator
          */
         virtual void setOperator(const EPMFloat factor, const Matrix& mat);

         /**
          * @brief Add operator's diagonal value
          *
          * @param factor Multiplicative factor
          */
         virtual void addOperator(const EPMFloat factor);

         /**
          * @brief Add operator value
          *
          * @param factor Multiplicative factor
          * @param mat    Input operator
          */
         virtual void addOperator(const EPMFloat factor, const Matrix& mat);

         /**
          * @brief Convert the stored operator into a bounded one
          */
         virtual void constructBoundedOperator();

         /**
          * @brief Solve linear equation with zero BC value
          *
          * @param vector RHS of the linear equation
          */
         virtual void solveZero(Array& vector);

      protected:
         Matrix   mBCRows;

      private:
   };

   template <typename TOpType> BCRowOperator<TOpType>::BCRowOperator(const int nBC, const int nTau, const int id)
      : BoundedOperatorBase<TOpType>(nBC, nTau, id), mBCRows(nBC, nTau)
   {
   }

   template <typename TOpType> inline void BCRowOperator<TOpType>::implementBCs(const Matrix& bcRows)
   {
      this->mBCRows = bcRows;
   }

   template <typename TOpType> inline void BCRowOperator<TOpType>::setOperator(const EPMFloat factor)
   {
      this->rOp().diagonal().setConstant(factor);
   }

   template <typename TOpType> inline void BCRowOperator<TOpType>::setOperator(const EPMFloat factor, const Matrix& mat)
   {
      this->rOp() = factor * mat;
   }

   template <typename TOpType> inline void BCRowOperator<TOpType>::addOperator(const EPMFloat factor)
   {
      this->rOp().diagonal().array() += factor;
   }

   template <typename TOpType> inline void BCRowOperator<TOpType>::addOperator(const EPMFloat factor, const Matrix& mat)
   {
      this->rOp() += factor * mat;
   }

   template <typename TOpType> void BCRowOperator<TOpType>::constructBoundedOperator()
   {
      // Impose the boundary condition on the last rows
      int start = this->nTau() - this->nBC();

      // Loop over all the boundary conditions
      for(int i = 0; i < this->nBC(); ++i)
      {
         this->rOp().row(start + i) = this->mBCRows.row(i);
      }
   }

   template <typename TOpType> void BCRowOperator<TOpType>::solveZero(Array& vector)
   {
      // Imposed on the boundary condition on the last rows
      vector.bottomRows(this->nBC()).setConstant(0.0);
      
      // Call basic solve
      this->solveEquation(vector);
   }

}

#endif // BCROWOPERATOR_HPP
