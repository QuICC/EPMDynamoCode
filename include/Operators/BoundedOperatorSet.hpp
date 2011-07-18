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

namespace EPMDynamo {

   /**
    * \brief Implementation of a set of bounded operators
    *
    * \tparam TOpType Type of the operator
    * \tparam TBOpType Type of the bounded operator
    *
    * \epmBug Protected part has beend made public. Is there a better solution?
    */
   template <typename TOpType, template <typename> class TBOpType> class BoundedOperatorSet
   {
      public:
         /// typedef for a homogeneous bounded operator
         typedef  TBOpType<TOpType> BOperator;

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

         /**
          * @brief Solve single equation with degree l operator with zero BC value
          *
          * @param rhs RHS and solution of equation
          * @param l Harmonic degree l
          */
         void solveVector(Array&  rhs, const int l, const bool isReal);
         
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
   };

   template <typename TOpType, template <typename> class TBOpType> BoundedOperatorSet<TOpType, TBOpType>::BoundedOperatorSet(SmartTruncation pTrunc, bool hasL0)
      : mHasL0(hasL0), mMinL(0), mNOp(0), mNbc(0), mpTrunc(pTrunc), mTmp()
   {
      // Set the minimum l index 
      if(this->mpTrunc->local()->spec()->lArray()(0) == 0)
      {
         this->mMinL = (! this->mHasL0);
      }
   }

   template <typename TOpType, template <typename> class TBOpType> inline bool BoundedOperatorSet<TOpType, TBOpType>::hasL0() const
   {
      return this->mHasL0;
   }

   template <typename TOpType, template <typename> class TBOpType> inline int BoundedOperatorSet<TOpType, TBOpType>::minL() const
   {
      return this->mMinL;
   }

   template <typename TOpType, template <typename> class TBOpType> inline int BoundedOperatorSet<TOpType, TBOpType>::nOp() const
   {
      return this->mNOp;
   }

   template <typename TOpType, template <typename> class TBOpType> inline int BoundedOperatorSet<TOpType, TBOpType>::nBC() const
   {
      return this->mNbc;
   }

   template <typename TOpType, template <typename> class TBOpType> inline SmartTruncation BoundedOperatorSet<TOpType, TBOpType>::trunc() const
   {
      return this->mpTrunc;
   }

   template <typename TOpType, template <typename> class TBOpType> inline void BoundedOperatorSet<TOpType, TBOpType>::addBC(SmartBC pBC) 
   {
      // Add BC
      this->mBCs.push_back(pBC);

      // Increment BC counter
      ++this->mNbc;
   }

   template <typename TOpType, template <typename> class TBOpType> inline const typename BoundedOperatorSet<TOpType, TBOpType>::BOperator& BoundedOperatorSet<TOpType, TBOpType>::harmOp(const int l) const
   {
      // Protect against use of unrequired (and probably dangerous) L=0 mode
      assert(!(this->mpTrunc->local()->spec()->lArray()(l) == 0 && (! this->mHasL0)));

      return *(this->mpOperators.at(l));
   }

   template <typename TOpType, template <typename> class TBOpType> inline typename BoundedOperatorSet<TOpType, TBOpType>::BOperator& BoundedOperatorSet<TOpType, TBOpType>::rHarmOp(const int l)
   {
      // Protect against use of unrequired (and probably dangerous) L=0 mode
      assert(!(this->mpTrunc->local()->spec()->lArray()(l) == 0 && (! this->mHasL0)));

      return *(this->mpOperators.at(l));
   }

   template <typename TOpType, template <typename> class TBOpType> void BoundedOperatorSet<TOpType, TBOpType>::initOperators()
   {
      // Get operator size
      int opSize = this->trunc()->sim()->rad()->nN();

       // Create BC rows matrix
      Matrix   bcRows(this->nBC(), opSize);

       // Create BC values array
      ArrayZ   bcVals(this->nBC());

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

            bcVals(k) = this->mBCs.at(k)->getRHSBC(l,this->trunc()->local()->spec()->m(0,l));
         }

         // Create smart pointer
         pBOp = SmartBOperator(new BOperator(this->nBC(), opSize, l));

         // Add homogeneous operator
         this->mpOperators.push_back(pBOp);

         // implement boundary condition of newly created operator
         this->mpOperators.back()->implementBCs(bcRows, bcVals);
      }

      // Set the total number of operators
      this->mNOp = this->mpOperators.size();
   }

   template <typename TOpType, template <typename> class TBOpType> void BoundedOperatorSet<TOpType, TBOpType>::multiplyOrders(MatrixZ& rhs, const MatrixZ& old, const int l)
   {
      // Get size of the matrix
      int rows = this->harmOp(l).nTau();
      int nM = this->trunc()->local()->spec()->nM(l);

      rhs.topLeftCorner(rows, nM) = this->harmOp(l).op()*old.topLeftCorner(rows,nM);
   }

   template <typename TOpType, template <typename> class TBOpType> void BoundedOperatorSet<TOpType, TBOpType>::affineOrders(MatrixZ& rhs, const MatrixZ& old, const MatrixZ& cTerms, const int l)
   {
      // Get size of the matrix
      int rows = this->harmOp(l).nTau();
      int nM = this->trunc()->local()->spec()->nM(l);

      rhs.topLeftCorner(rows, nM) = this->harmOp(l).op()*old.topLeftCorner(rows,nM) + cTerms.topLeftCorner(rows, nM);
   }

   template <typename TOpType, template <typename> class TBOpType> inline void BoundedOperatorSet<TOpType, TBOpType>::solveZeroOrders(MatrixZ& rhs, const int l)
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

   template <typename TOpType, template <typename> class TBOpType> inline void BoundedOperatorSet<TOpType, TBOpType>::solveOrders(MatrixZ& rhs, const int l)
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
         this->rHarmOp(l).solve(this->mTmp, true);

         // Copy solution into field
         for(int j = 0; j < nVar; ++j)
         {
            rhs(j,m).real() = this->mTmp(j);
         }

         // Set tmp variable to imaginary part
         this->mTmp = rhs.col(m).imag();

         // Solve imaginary equation
         this->rHarmOp(l).solve(this->mTmp, false);

         // Copy solution into field
         for(int j = 0; j < nVar; ++j)
         {
            rhs(j,m).imag() = this->mTmp(j);
         }
      }
   }

   template <typename TOpType, template <typename> class TBOpType> inline void BoundedOperatorSet<TOpType, TBOpType>::solveZeroVector(Array& rhs, const int l)
   {
      this->rHarmOp(l).solveZero(rhs);
   }

   template <typename TOpType, template <typename> class TBOpType> inline void BoundedOperatorSet<TOpType, TBOpType>::solveVector(Array& rhs, const int l, const bool isReal)
   {
      this->rHarmOp(l).solve(rhs, isReal);
   }
}

#endif // BOUNDEDOPERATORSET_HPP
