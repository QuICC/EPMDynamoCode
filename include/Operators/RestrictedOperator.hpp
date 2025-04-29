/** \file RestrictedOperator.hpp
 *  \brief Implementation of a restricted bounded operator
 */

#ifndef RESTRICTEDOPERATOR_HPP
#define RESTRICTEDOPERATOR_HPP

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
    * \brief Implementation of a restricted bounded operator
    *
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class RestrictedOperator: public BoundedOperatorBase<TOpType>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nBC Number of boundary conditions
          * @param nBOp Size of the full operator
          * @param id ID of the operator (for example harmonic degree l)
          * @param simID Simulation wide ID of the operator (for example harmonic degree l)
          */
         RestrictedOperator(const int nBC, const int nNOp, const int id, const int simID);

         /**
          * @brief Destructor
          */
         virtual ~RestrictedOperator() {};

         /**
          * @brief Implement the boundary conditions
          */
         virtual void implementBCs();

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

         /**
          * @brief Solve linear equation
          *
          * @param vector RHS of the linear equation
          *
          * \epmBug THIS IMPLEMENTATION IS WRONG
          */
         virtual void solve(Array& vector, const int simM, const bool isReal);

         /**
          * @brief Extend solution to full truncation
          */
         void extendZero(MatrixZ& rMat) const;

         /**
          * @brief Extend solution to full truncation
          */
         void extendZero(Array& rKernel) const;

      protected:

      private:
         /**
          * @brief Extension matrix storage
          */
         ExtensionMatrix  mExt;
   };

   template <typename TOpType> RestrictedOperator<TOpType>::RestrictedOperator(const int nBC, const int nNOp, const int id, const int simID)
      : BoundedOperatorBase<TOpType>(nBC, nNOp-nBC, id, simID), mExt()
   {
   }

   template <typename TOpType> inline void RestrictedOperator<TOpType>::implementBCs()
   {
       // Create BC rows matrix
      Matrix   bcRows(this->nBC(), this->nTau()+this->nBc());

      // Create boundary rows from boundary conditions
      for(int k = 0; k < this->nBC(); ++k)
      {
         bcRows.row(k) = this->mBCs.at(k)->getLHSBC(this->id()).transpose();
      }

      // Build the extension matrices
      this->mExt.build(bcRows);
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::setOperator(const EPMFloat factor)
   {
      this->mExt.setRestrictedOperator(this->rOp(), factor, Matrix::Identity(this->nTau()+this->nBC(),this->nTau()+this->nBC()));
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::setOperator(const EPMFloat factor, const Matrix& mat)
   {
      this->mExt.setRestrictedOperator(this->rOp(), factor, mat);
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::addOperator(const EPMFloat factor)
   {
      this->mExt.addRestrictedOperator(this->rOp(), factor, Matrix::Identity(this->nTau()+this->nBC(),this->nTau()+this->nBC()));
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::addOperator(const EPMFloat factor, const Matrix& mat)
   {
      this->mExt.addRestrictedOperator(this->rOp(), factor, mat);
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::constructBoundedOperator()
   {
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::solveZero(Array& vector)
   {
      // Call basic solve
      this->solveEquation(vector);

      // Extend solution to full truncation
      this->mExt.extendZero(vector);
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::solve(Array& vector, const int simM, const bool isReal)
   {
      // Call basic solve
      this->solveEquation(vector);

      // Extend solution to full truncation
      this->mExt.extendZero(vector);
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::extendZero(MatrixZ& rMat) const
   {
      // Extend solution to full truncation
      this->mExt.extendZero(rMat);
   }

   template <typename TOpType> void RestrictedOperator<TOpType>::extendZero(Array& rKernel) const
   {
      // Extend solution to full truncation
      this->mExt.extendZero(rKernel);
   }

}

#endif // RESTRICTEDOPERATOR_HPP
