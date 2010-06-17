/** \file HomogeneousBOperator.hpp
 *  \brief Implementation of a bounded operator on a set of homogeneous boundary conditions
 */

#ifndef HOMOGENEOUSBOPERATOR_HPP
#define HOMOGENEOUSBOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Operators/BoundedOperator.hpp"
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/PolyOps/POperator.hpp"
#include "LAPACK_Iface.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a bounded operator on a set of homogeneous boundary conditions
    *
    * @param TOpType Type of the operator
    *
    * \bug Extension matrix computation has to be restructured/corrected
    */
   template <typename TOpType> class HomogeneousBOperator: public BoundedOperator<TOpType>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nBC Number of boundary conditions
          * @param size Size of the operator (size of the basis)
          * @param id ID of the operator (for example harmonic degree)
          */
         HomogeneousBOperator(const int nBC, const int size, const int id);

         /**
          * @brief Destructor
          */
         virtual ~HomogeneousBOperator() {};

         /**
          * @brief Implement the boundary conditions
          *
          * @param bcRows Matrix of the boundary imposing rows
          */
         void implementBCs(const Matrix& bcRows);

         /**
          * @brief Compute the operator the bounded operator based on input
          *
          * @param factor Multiplicative factor
          * @param polyOp Polynomial operator defining the operator
          *
          * \bug PolynomialOperators have to be restructured
          */
         void constructBOperator(const DynamoFloat factor, const PolynomialOperator<POperator>& polyOp);

         /**
          * @brief Solve for vector
          *
          * @param vector RHS of the linear equation
          */
         void solve(Array& vector);
         
      protected:
         /**
          * @brief Get the extension matrix
          */
         const SparseMatrix&  extE() const;

         /**
          * @brief Set the extension matrix
          */
         SparseMatrix&  rExtE();

         /**
          * @brief Compute extension matrix
          *
          * @param bcRows BC imposing rows
          */
         void createExtension(Matrix& bcRows);

      private:
         /**
          * @brief Extension matrix storage
          */
         SparseMatrix   mE;

         /**
          * @brief Convert BCRows into extension matrix rows
          *
          * @param bcRows Convert BC rows to extension matrix rows
          */
         void convertBCRows(Matrix& bcRows);
   };

   template <typename TOpType> inline const SparseMatrix& HomogeneousBOperator<TOpType>::extE() const
   {
      return this->mE;
   }

   template <typename TOpType> inline SparseMatrix& HomogeneousBOperator<TOpType>::rExtE()
   {
      return this->mE;
   }

   template <typename TOpType> HomogeneousBOperator<TOpType>::HomogeneousBOperator(const int nBC, const int size, const int id)
      : BoundedOperator<TOpType>(nBC, size-nBC, id), mE(size, size-nBC)
   {
   }

   template <typename TOpType> inline void HomogeneousBOperator<TOpType>::implementBCs(const Matrix& bcRows)
   {
      // Copy boundary rows to temporary
      Matrix tmp = bcRows;

      // Compute extension matrix
      this->createExtension(tmp);
   }

   template <typename TOpType> void HomogeneousBOperator<TOpType>::constructBOperator(const DynamoFloat factor, const PolynomialOperator<POperator>& polyOp)
   {
      // Get sizes
      int rows = this->nTau();
      int cols = this->extE().rows();

      // Compute bounded operator
      this->rOp() = factor*polyOp.op().corner(Eigen::TopLeft, rows, cols) * this->extE();
   }

   template <typename TOpType> void HomogeneousBOperator<TOpType>::solve(Array& vector)
   {
      // Call basic solve
      this->solveEquation(vector);

      // fill in value for last coefficient(s)
      vector = this->extE()*vector.start(this->nTau());
   }

   template <typename TOpType> void HomogeneousBOperator<TOpType>::createExtension(Matrix& bcRows)
   {
      // Compute the BC imposing part of the matrix
      this->convertBCRows(bcRows);

      // Get sizes
      int rows = this->extE().rows();
      int cols = this->nTau();

      // Initialise sparse extension matrix
      this->rExtE().startFill(cols + this->nBC()*rows); // estimate of the number of nonzeros

      // Fill sparse extension matrix
      for (int j=0; j < cols; ++j)
      {
         for (int i=0; i < rows; ++i)
         {
            // Fill identity part
            if(i < cols && i == j)
            {
               this->rExtE().fill(i,j) = 1.0;
            }
            // Fill BC part
            else if(i >= cols)
            {
               this->rExtE().fill(i,j) = bcRows(i-cols,j);
            }
         }
      }
      // End filling of sparse matrix
      this->rExtE().endFill();
   }

   template <typename TOpType> void HomogeneousBOperator<TOpType>::convertBCRows(Matrix& bcRows)
   {
      // Get matrix size
      const int nN = this->extE().rows();
      // Create temporary matrix and set to idendity
      Matrix tmp(nN, nN);
      tmp.setIdentity();

      // Get start row of BC values
      int bcStart = this->nTau();

      if(this->nBC() > 0)
      {
         // Set boundary rows
         for(int i = bcStart; i < nN; ++i)
         {
            tmp.row(i) = bcRows.row(i-bcStart);
         }

         // Compute LU factorisation
         int rows = tmp.rows();
         int cols = tmp.cols();
         int info;
         ArrayI ipiv(std::min(rows, cols));
         dgetrf_(&rows, &cols, tmp.data(), &rows, ipiv.data(), &info);
      
         // Compute Inverse of matrix
         int lwork = 5*rows;
         Array work(lwork);
         dgetri_(&rows, tmp.data(), &rows, ipiv.data(), work.data(), &lwork, &info);
      }

      // Put result back into bcRows
      for(int i=0; i < this->nBC(); ++i)
      {
         bcRows.row(i) = tmp.row(bcStart+i);
      }
   }

}

#endif // HOMOGENEOUSBOPERATOR_HPP
