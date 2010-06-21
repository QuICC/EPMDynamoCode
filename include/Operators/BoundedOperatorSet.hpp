/** \file BoundedOperatorSet.hpp
 *  \brief Implementation of a set of bounded operators
 */

#ifndef BOUNDEDOPERATORSET_HPP
#define BOUNDEDOPERATORSET_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//
#include <vector>
#include <set>
#include <map>

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Operators/BoundedOperatorBase.hpp"
#include "Operators/BoundedOperator.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a set of bounded operators
    *
    * \tparam TOpType Type of the operator
    *
    * \bug Need to sort out the redirect, Real/Imag, m values issues
    * \bug The operators are not yet computed correctly
    */
   template <typename TOpType> class BoundedOperatorSet
   {
      public:
         /// typedef for a homogeneous bounded operator
         typedef  BoundedOperator<TOpType, Homogeneous> HBOperator;

         /// typedef for a non homogeneous bounded operator
         typedef  BoundedOperator<TOpType, NonHomogeneous> NHBOperator;

         /// typedef for a base bounded operator
         typedef  BoundedOperatorBase<TOpType> BOperatorBase;

         /// typedef for a smart pointer to a base bounded operator
         typedef  EPMSHARED_PTR<BOperatorBase> SmartBOperatorBase;

         /// typedef for a smart pointer to a homogeneous bounded operator
         typedef  EPMSHARED_PTR<HBOperator> SmartHBOperator;

         /// typedef for a smart pointer to a non homogeneous bounded operator
         typedef  EPMSHARED_PTR<HBOperator> SmartNHBOperator;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         BoundedOperatorSet(SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~BoundedOperatorSet() {};

         /**
          * @brief Init the operators
          */
         virtual void initOperators();

         /**
          * @brief Add a BoundaryCondition pointer
          *
          * @param pBC Smart pointer to boundary condition
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Get truncation information
          */
         SmartTruncation trunc() const;

         /**
          * @brief Compute affine the matrix product for the orders of given degree
          *
          * @param rhs RHS output value
          * @param old Old RHS value
          * @param cTerms Non linear terms to add
          * @param l Harmonic degree l
          *
          * \bug Is this required like this?
          */
         void affineOrders(MatrixZ&  rhs, const MatrixZ& old, const MatrixZ& cTerms, const int l);

         /**
          * @brief Solve equations for all orders with degree l
          *
          * \bug THIS IS JUST A TEMPORARY WORKAROUNG, but to what?
          * \bug Is this required like this?
          */
         void solveVector(Array&  rhs, const int l);
         
      protected:
         /**
          * @brief Get the number of operators
          */
         int nOp() const;

         /**
          * @brief Get the number of boundary conditions
          */
         int nBC() const;

         /**
          * @brief Compute the matrix product for the orders of given degree
          *
          * @param rhs RHS output
          * @param old Old RHS values
          * @param l Harmonic degree l
          *
          * \bug Is this required like this?
          */
         void multiplyOrders(MatrixZ&  rhs, const MatrixZ& old, const int l);

         /**
          * @brief Solve equations for all orders with degree l
          *
          * @param Equation rhs and solution
          * @param l Harmonic degree l
          *
          * \bug Is this required like this?
          */
         void solveOrders(MatrixZ&  rhs, const int l);

         /**
          * @brief Solve equations for all orders with degree l with zero BC value
          *
          * @param rhs Equation rhs and solution
          * @param l Harmonic degree l
          *
          * \bug Is this required like this?
          */
         void solveZeroOrders(MatrixZ&  rhs, const int l);

         /**
          * @brief Vector of BoundaryCondition
          */
         std::vector<SmartBC>    mBCs;

         /**
          * @brief Get directe access to operator without redirection
          *
          * @param i Index of operator
          */
         const BOperatorBase& op(const int i) const;

         /**
          * @brief Set directe access to operator without redirection
          *
          * @param i Index of operator
          */
         BOperatorBase& rOp(const int i);

      private:
         /**
          * @brief Number of operators in the set
          */
         int mNOp;

         /**
          * @brief Number of boundary conditions
          */
         int mNbc;

         /**
          *  @brief Truncation information
          */
         SmartTruncation mpTrunc;

         /**
          * @brief An array for temporary storage
          */
         Array mTmp;

         /**
          * @brief Vector of operators
          */
         std::vector<SmartBOperatorBase> mpOperators;

         /**
          * @brief Convert l,m index to vector of operators index
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          *
          *\bug Redirection or the general aspect of minimal memory implementation has to be defined first. This is currently wrong
          */
         int redirect(const int l, const int m) const;

         /**
          * @brief Get Operator for real component for harmonic l,m
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          */
         const BOperatorBase& realOp(const int l, const int m) const;

         /**
          * @brief Get Operator for imaginary component for harmonic l,m
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          */
         const BOperatorBase& imagOp(const int l, const int m) const;

         /**
          * @brief Set Operator for real component for harmonic l,m
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          */
         BOperatorBase& rRealOp(const int l, const int m);

         /**
          * @brief Set Operator for imaginary component for harmonic l,m
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          */
         BOperatorBase& rImagOp(const int l, const int m);
   };

   template <typename TOpType> BoundedOperatorSet<TOpType>::BoundedOperatorSet(SmartTruncation pTrunc)
      : mNOp(0), mNbc(0), mpTrunc(pTrunc), mTmp()
   {
   }

   template <typename TOpType> inline int BoundedOperatorSet<TOpType>::nOp() const
   {
      return this->mNOp;
   }

   template <typename TOpType> inline int BoundedOperatorSet<TOpType>::nBC() const
   {
      return this->mNbc;
   }

   template <typename TOpType> inline SmartTruncation BoundedOperatorSet<TOpType>::trunc() const
   {
      return this->mpTrunc;
   }

   template <typename TOpType> inline void BoundedOperatorSet<TOpType>::addBC(SmartBC pBC) 
   {
      // Add BC
      this->mBCs.push_back(pBC);

      // Increment BC counter
      ++this->mNbc;
   }

   template <typename TOpType> inline int BoundedOperatorSet<TOpType>::redirect(const int l, const int m)  const
   {
      return l;
   }

   template <typename TOpType> inline const typename BoundedOperatorSet<TOpType>::BOperatorBase& BoundedOperatorSet<TOpType>::realOp(const int l, const int m) const
   {
      return *(this->mpOperators.at(this->redirect(l,m)));
   }

   template <typename TOpType> inline const typename BoundedOperatorSet<TOpType>::BOperatorBase& BoundedOperatorSet<TOpType>::imagOp(const int l, const int m) const
   {
      return *(this->mpOperators.at(this->redirect(l,m)));
   }

   template <typename TOpType> inline typename BoundedOperatorSet<TOpType>::BOperatorBase& BoundedOperatorSet<TOpType>::rRealOp(const int l, const int m)
   {
      return *(this->mpOperators.at(this->redirect(l,m)));
   }

   template <typename TOpType> inline typename BoundedOperatorSet<TOpType>::BOperatorBase& BoundedOperatorSet<TOpType>::rImagOp(const int l, const int m)
   {
      return *(this->mpOperators.at(this->redirect(l,m)));
   }

   template <typename TOpType> inline const typename BoundedOperatorSet<TOpType>::BOperatorBase& BoundedOperatorSet<TOpType>::op(const int i) const
   {
      return *(this->mpOperators.at(i));
   }

   template <typename TOpType> inline typename BoundedOperatorSet<TOpType>::BOperatorBase& BoundedOperatorSet<TOpType>::rOp(const int i)
   {
      return *(this->mpOperators.at(i));
   }

   template <typename TOpType> void BoundedOperatorSet<TOpType>::initOperators()
   {
      // Get operator size
      int opSize = this->trunc()->sim()->rad()->nN();
       // Create BC rows matrix
      Matrix   bcRows(this->nBC(), opSize);

      // Get number of harmonic degrees
      int nL = this->trunc()->local()->spec()->nL();

      // Pointers to homogeneous operators
      SmartHBOperator   pHOp;

      int nM;

      // Loop over degreees
      for(int l = 0; l < nL; ++l)
      {
         nM = this->trunc()->local()->spec()->nM(l);
         for(int m = 0; m < nM; ++m)
         {
            assert(this->mBCs.at(0)->getRealRHSBC(l,m) == 0.0);
            assert(this->mBCs.at(0)->getImagRHSBC(l,m) == 0.0);
         }

         // Get boundary row
         bcRows.row(0) = this->mBCs.at(0)->getLHSBC(l).transpose();

         // Create smart pointer
         pHOp = SmartHBOperator(new HBOperator(this->nBC(), opSize, l));

         // Add homogeneous operator
         this->mpOperators.push_back(pHOp);

         // implement boundary condition of newly created operator
         this->mpOperators.back()->implementBCs(bcRows);
      }

      // Set the total number of operators
      this->mNOp = this->mpOperators.size();
   }

   template <typename TOpType> void BoundedOperatorSet<TOpType>::multiplyOrders(MatrixZ& rhs, const MatrixZ& old, const int l)
   {
      // Get size of the matrix
      int rows;
      int nM = this->trunc()->local()->spec()->nM(l);

      // Loop over the orders
      for(int m = 0; m < nM; ++m)
      {
         // Get size of real operator
         rows = this->realOp(l,m).nTau();

         // Put real values in tmp
         this->mTmp = this->realOp(l,m).op() * old.col(m).start(rows).real();

         // Copy tmp into real part of rhs
         for(int j=0; j < rows; ++j)
         {
            rhs(j,m).real() = this->mTmp(j);
         }

         // Get size of imaginary operator
         rows = this->imagOp(l,m).nTau();

         // Put imaginary values in tmp
         this->mTmp = this->imagOp(l,m).op() * old.col(m).start(rows).imag();

         // Copy tmp into imaginary part of rhs
         for(int j=0; j < rows; ++j)
         {
            rhs(j,m).imag() = this->mTmp(j);
         }
      }
   }

   template <typename TOpType> void BoundedOperatorSet<TOpType>::affineOrders(MatrixZ& rhs, const MatrixZ& old, const MatrixZ& cTerms, const int l)
   {
      // Get size of the matrix
      int rows;
      int nM = this->trunc()->local()->spec()->nM(l);

      // Loop over the orders
      for(int m = 0; m < nM; ++m)
      {
         // Get size of real operator
         rows = this->realOp(l,m).nTau();

         // Put real values in tmp
         this->mTmp = this->realOp(l,m).op() * old.col(m).start(rows).real() + cTerms.col(m).start(rows).real();

         // Copy tmp into real part of rhs
         for(int j=0; j < rows; ++j)
         {
            rhs(j,m).real() = this->mTmp(j);
         }

         // Get size of real operator
         rows = this->imagOp(l,m).nTau();

         // Put imaginary values in tmp
         this->mTmp = this->imagOp(l,m).op() * old.col(m).start(rows).imag() + cTerms.col(m).start(rows).imag();

         // Copy tmp into imaginary part of rhs
         for(int j=0; j < rows; ++j)
         {
            rhs(j,m).imag() = mTmp(j);
         }
      }
   }

   template <typename TOpType> void BoundedOperatorSet<TOpType>::solveOrders(MatrixZ& rhs, const int l)
   {
      // Get size of the matrix
      int nVar = rhs.rows();
      int nM = this->trunc()->local()->spec()->nM(l);

      // Loop over the orders
      for(int m = 0; m < nM; ++m)
      {
         // Set tmp variable to real part
         this->mTmp = rhs.col(m).real();

         // Solve real equation
         this->rRealOp(l,m).solve(this->mTmp);

         // Copy solution into field
         for(int j = 0; j < nVar; ++j)
         {
            rhs(j,m).real() = this->mTmp(j);
         }

         // Set tmp variable to imaginary part
         this->mTmp = rhs.col(m).imag();

         // Solve imaginary equation
         this->rImagOp(l,m).solve(this->mTmp);

         // Copy solution into field
         for(int j = 0; j < nVar; ++j)
         {
            rhs(j,m).imag() = this->mTmp(j);
         }
      }
   }

   template <typename TOpType> inline void BoundedOperatorSet<TOpType>::solveZeroOrders(MatrixZ& rhs, const int l)
   {
afafafafaf;
   }

   template <typename TOpType> inline void BoundedOperatorSet<TOpType>::solveVector(Array& rhs, const int l)
   {
      this->rRealOp(l,0).solve(rhs);
   }

}

#endif // BOUNDEDOPERATORSET_HPP
