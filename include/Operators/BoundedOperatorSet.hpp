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
    * \epmBug Protected part has beend made public. Is there a better solution?
    */
   template <typename TOpType> class BoundedOperatorSet
   {
      public:
         /// typedef for a homogeneous bounded operator
         typedef  BoundedOperator<TOpType> BOperator;

         /// typedef for a smart pointer to a bounded operator
         typedef  EPMSHARED_PTR<BOperator> SmartBOperator;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param hasL0 Is the L = 0 mode required ?
          */
         BoundedOperatorSet(SmartTruncation pTrunc, bool hasL0);

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
          */
         void affineOrders(MatrixZ&  rhs, const MatrixZ& old, const MatrixZ& cTerms, const int l);

         /**
          * @brief Solve single equation with degree l operator with zero BC value
          *
          * @param rhs RHS and solution of equation
          * @param l Harmonic degree l
          */
         void solveZeroVector(Array&  rhs, const int l);
         
//      protected:
         /**
          * @brief Is the l=0 mode required?
          */
         bool hasL0() const;

         /**
          * @brief Minimul harmonic degree index
          */
         int minL() const;

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
          */
         void multiplyOrders(MatrixZ&  rhs, const MatrixZ& old, const int l);

         /**
          * @brief Solve equations for all orders with degree l
          *
          * @param rhs Equation rhs and solution
          * @param l Harmonic degree l
          */
         void solveOrders(MatrixZ&  rhs, const int l);

         /**
          * @brief Solve equations for all orders with degree l with zero BC value
          *
          * @param rhs Equation rhs and solution
          * @param l Harmonic degree l
          */
         void solveZeroOrders(MatrixZ&  rhs, const int l);

         /**
          * @brief Vector of BoundaryCondition
          */
         std::vector<SmartBC>    mBCs;

         /**
          * @brief Get operator for given harmonic degree index
          *
          * @param l Harmonic degree index
          */
         const BOperator& harmOp(const int l) const;

         /**
          * @brief Set operator for given harmonic degree index
          *
          * @param l Harmonic degree index
          */
         BOperator& rHarmOp(const int l);

      private:
         /**
          * @brief Is the l=0 mode required ?
          */
         bool mHasL0;

         /**
          * @brief Minimum harmonic degree index (depends on mHasL0)
          */
         int mMinL;

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
         std::vector<SmartBOperator> mpOperators;

         /**
          * @brief Flag to check if all BCs are homogeneous
          */
         bool mAllHomogeneous;

         /**
          * @brief Check if the NH matrix is required for harmonic degree l
          *
          * @param l Harmonic degree l
          */
         bool needNHMatrix(const int l) const;

         /**
          * @brief Check if real boundary conditions are non homogeneous
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          */
         bool isRealBCNH(const int l, const int m) const;

         /**
          * @brief Check if real boundary conditions are non homogeneous
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          */
         bool isImagBCNH(const int l, const int m) const;
   };

   template <typename TOpType> BoundedOperatorSet<TOpType>::BoundedOperatorSet(SmartTruncation pTrunc, bool hasL0)
      : mHasL0(hasL0), mMinL(0), mNOp(0), mNbc(0), mpTrunc(pTrunc), mTmp(), mAllHomogeneous(false)
   {
      // Set the minimum l index 
      if(this->mpTrunc->local()->spec()->lArray()(0) == 0)
      {
         this->mMinL = (! this->mHasL0);
      }
   }

   template <typename TOpType> inline bool BoundedOperatorSet<TOpType>::hasL0() const
   {
      return this->mHasL0;
   }

   template <typename TOpType> inline int BoundedOperatorSet<TOpType>::minL() const
   {
      return this->mMinL;
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

   template <typename TOpType> inline const typename BoundedOperatorSet<TOpType>::BOperator& BoundedOperatorSet<TOpType>::harmOp(const int l) const
   {
      // Protect against use of unrequired (and probably dangerous) L=0 mode
      assert(!(this->mpTrunc->local()->spec()->lArray()(l) == 0 && (! this->mHasL0)));

      return *(this->mpOperators.at(l));
   }

   template <typename TOpType> inline typename BoundedOperatorSet<TOpType>::BOperator& BoundedOperatorSet<TOpType>::rHarmOp(const int l)
   {
      // Protect against use of unrequired (and probably dangerous) L=0 mode
      assert(!(this->mpTrunc->local()->spec()->lArray()(l) == 0 && (! this->mHasL0)));

      return *(this->mpOperators.at(l));
   }

   template <typename TOpType> void BoundedOperatorSet<TOpType>::initOperators()
   {
      // Get operator size
      int opSize = this->trunc()->sim()->rad()->nN();
      int BOpSize = opSize - this->nBC();
       // Create BC rows matrix
      Matrix   bcRows(this->nBC(), opSize);

      // Get number of harmonic degrees
      int nL = this->trunc()->local()->spec()->nL();

      // Pointers to homogeneous operators
      SmartBOperator   pBOp;

      // Loop over degreees
      for(int l = 0; l < nL; ++l)
      {
         // Get boundary rows
         for(int k=0; k < this->nBC(); ++k)
         {
            bcRows.row(k) = this->mBCs.at(k)->getLHSBC(l).transpose();
         }

         // Create smart pointer
         pBOp = SmartBOperator(new BOperator(this->nBC(), BOpSize, l));

         // Add homogeneous operator
         this->mpOperators.push_back(pBOp);

         // implement boundary condition of newly created operator
         this->mpOperators.back()->implementBCs(bcRows, this->needNHMatrix(l));
      }

      // Set the total number of operators
      this->mNOp = this->mpOperators.size();

      // Check if all BCs are homogeneous
      this->mAllHomogeneous = true;
      for(int k=0; k < this->nBC(); ++k)
      {
         this->mAllHomogeneous = this->mAllHomogeneous && this->mBCs.at(k)->isHomogeneous();
      }
   }

   template <typename TOpType> void BoundedOperatorSet<TOpType>::multiplyOrders(MatrixZ& rhs, const MatrixZ& old, const int l)
   {
      // Get size of the matrix
      int rows = this->harmOp(l).nTau();
      int nM = this->trunc()->local()->spec()->nM(l);

      // Loop over the orders
      for(int m = 0; m < nM; ++m)
      {
         // Put real values in tmp
         this->mTmp = this->harmOp(l).op() * old.col(m).start(rows).real();

         // Copy tmp into real part of rhs
         for(int j=0; j < rows; ++j)
         {
            rhs(j,m).real() = this->mTmp(j);
         }

         // Put imaginary values in tmp
         this->mTmp = this->harmOp(l).op() * old.col(m).start(rows).imag();

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
      int rows = this->harmOp(l).nTau();
      int nM = this->trunc()->local()->spec()->nM(l);

      // Loop over the orders
      for(int m = 0; m < nM; ++m)
      {
         // Put real values in tmp
         this->mTmp = this->harmOp(l).op() * old.col(m).start(rows).real() + cTerms.col(m).start(rows).real();

         // Copy tmp into real part of rhs
         for(int j=0; j < rows; ++j)
         {
            rhs(j,m).real() = this->mTmp(j);
         }

         // Put imaginary values in tmp
         this->mTmp = this->harmOp(l).op() * old.col(m).start(rows).imag() + cTerms.col(m).start(rows).imag();

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

      if(this->mAllHomogeneous)
      {
         // Loop over the orders
         for(int m = 0; m < nM; ++m)
         {
            // Set tmp variable to real part
            this->mTmp = rhs.col(m).real();

            // Solve real equation
            this->rHarmOp(l).solve(this->mTmp, true);

            // Copy solution into field
            for(int j = 0; j < nVar; ++j)
            {
               rhs(j,m).real() = this->mTmp(j);
            }

            // Set tmp variable to imaginary part
            this->mTmp = rhs.col(m).imag();

            // Solve imaginary equation
            this->rHarmOp(l).solve(this->mTmp, true);

            // Copy solution into field
            for(int j = 0; j < nVar; ++j)
            {
               rhs(j,m).imag() = this->mTmp(j);
            }
         }
      } else
      {
         // Loop over the orders
         for(int m = 0; m < nM; ++m)
         {
            // Set tmp variable to real part
            this->mTmp = rhs.col(m).real();

            // Include the RHS boundary value
            if(this->isRealBCNH(l, m))
            {
               for(int k = 0; k < this->nBC(); ++k)
               {
                  this->mTmp(nVar-this->nBC()+k) = this->mBCs.at(k)->getRealRHSBC(l,m);
               }
            }

            // Solve real equation
            this->rHarmOp(l).solve(this->mTmp, true);

            // Copy solution into field
            for(int j = 0; j < nVar; ++j)
            {
               rhs(j,m).real() = this->mTmp(j);
            }

            // Set tmp variable to imaginary part
            this->mTmp = rhs.col(m).imag();

            // Include the RHS boundary value
            if(this->isImagBCNH(l, m))
            {
               for(int k = 0; k < this->nBC(); ++k)
               {
                  this->mTmp(nVar-this->nBC()+k) = this->mBCs.at(k)->getImagRHSBC(l,m);
               }
            }

            // Solve imaginary equation
            this->rHarmOp(l).solve(this->mTmp, true);

            // Copy solution into field
            for(int j = 0; j < nVar; ++j)
            {
               rhs(j,m).imag() = this->mTmp(j);
            }
         }
      }
   }

   template <typename TOpType> inline void BoundedOperatorSet<TOpType>::solveZeroOrders(MatrixZ& rhs, const int l)
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
         this->rHarmOp(l).solveZero(this->mTmp);

         // Copy solution into field
         for(int j = 0; j < nVar; ++j)
         {
            rhs(j,m).real() = this->mTmp(j);
         }

         // Set tmp variable to imaginary part
         this->mTmp = rhs.col(m).imag();

         // Solve imaginary equation
         this->rHarmOp(l).solveZero(this->mTmp);

         // Copy solution into field
         for(int j = 0; j < nVar; ++j)
         {
            rhs(j,m).imag() = this->mTmp(j);
         }
      }
   }

   template <typename TOpType> inline void BoundedOperatorSet<TOpType>::solveZeroVector(Array& rhs, const int l)
   {
      this->rHarmOp(l).solveZero(rhs);
   }

   template <typename TOpType> bool BoundedOperatorSet<TOpType>::isRealBCNH(const int l, const int m) const
   {
      bool flag = true;

      for(int k = 0; k < this->nBC(); ++k)
      {
         flag = flag && (this->mBCs.at(k)->getRealRHSBC(l,m) == 0);
      }

      return flag;
   }

   template <typename TOpType> bool BoundedOperatorSet<TOpType>::isImagBCNH(const int l, const int m) const
   {
      bool flag = true;

      for(int k = 0; k < this->nBC(); ++k)
      {
         flag = flag && (this->mBCs.at(k)->getImagRHSBC(l,m) == 0);
      }

      return flag;
   }

   template <typename TOpType> bool BoundedOperatorSet<TOpType>::needNHMatrix(const int l) const
   {
      int nM = this->trunc()->local()->spec()->nM(l);
      bool flag = true;

      for(int m=0; m < nM; ++m)
      {
         flag = flag && (! this->isRealBCNH(l,m)) && (! this->isImagBCNH(l,m));
      }

      return flag;
   }
}

#endif // BOUNDEDOPERATORSET_HPP
