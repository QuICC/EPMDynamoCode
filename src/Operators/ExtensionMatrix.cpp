/** \file ExtensionMatrix.cpp
 *  \brief Source of the implementation of the extension matrix to impose boundary conditions
 */

// System includes
//

// External includes
//

// Class include
//
#include "Operators/ExtensionMatrix.hpp"

// Project includes
//

namespace EPMDynamo {

   ExtensionMatrix::ExtensionMatrix()
      : mNBCs(0), mNBOp(0)
   {
   }

   void ExtensionMatrix::build(const Matrix &bcRows)
   {
      // Get number of boundary conditions
      this->mNBCs = bcRows.rows();

      // Get size of bounded operator
      this->mNBOp = bcRows.cols() - this->mNBCs;

      // This resize the matrix G
      this->mGMatrix.resize(this->mNBCs, this->mNBOp);

      // This resize the matrix H
      this->mHMatrix.resize(this->mNBCs, this->mNBCs);

      // Fill in the values for H
      this->mHMatrix = bcRows.bottomRightCorner(this->mNBCs, this->mNBCs);

      // Fill in the values for G
      this->mGMatrix = bcRows.topLeftCorner(this->mNBCs, this->mNBOp);

      // Invert the content of H to get the real H matrix
      this->invertHMatrix();
      
      // Compute matrix G
      this->mGMatrix = -this->H() * this->G();
   }

   void ExtensionMatrix::invertHMatrix()
   {
      // Get LU factorisation of matrix
      int rows = this->mNBCs;
      int cols = this->mNBCs;

      int info;
      ArrayI ipiv(rows);

      // Call LAPACK dgetrf routine for factorisation
      dgetrf_(&rows, &cols, this->mHMatrix.data(), &rows, ipiv.data(), &info);

      // Test success of computation through assert
      assert(info == 0);

      // Prepare for inversion
      rows = this->mNBCs;
      int lwork = -1;
      Array work(rows);

      // Call LAPACK dgetri routine for optimal work size
      dgetri_(&rows, this->mHMatrix.data(), &rows, ipiv.data(), work.data(), &lwork, &info);

      // Test success of computation through assert
      assert(info == 0);

      // Set to optimal work size
      lwork = static_cast<int>(work(0));
      work.resize(lwork);

      // Call LAPACK dgetri routine for inverse computation
      dgetri_(&rows, this->mHMatrix.data(), &rows, ipiv.data(), work.data(), &lwork, &info);

      // Test success of computation through assert
      assert(info == 0);
   }

   void ExtensionMatrix::prepareRHS(Array &rVector)
   {
      rVector.head(this->mNBOp) -= this->NH() * rVector.tail(this->mNBCs);
   }

   void ExtensionMatrix::extendZero(Array &rVector) const
   {
      rVector.tail(this->mNBCs) = this->G() * rVector.head(this->mNBOp);
   }

   void ExtensionMatrix::extendZero(MatrixZ &rMat) const
   {
      rMat.bottomRows(this->mNBCs) = this->G() * rMat.topRows(this->mNBOp);
   }

   void ExtensionMatrix::extend(Array &rVector) const
   {
      rVector.tail(this->mNBCs) = this->G() * rVector.head(this->mNBOp) + this->H() * rVector.tail(this->mNBCs);
   }

   void ExtensionMatrix::extend(MatrixZ &rMat) const
   {
      rMat.bottomRows(this->mNBCs) = this->G() * rMat.topRows(this->mNBOp) + this->H() * rMat.bottomRows(this->mNBCs);
   }

   void ExtensionMatrix::setRestrictedOperator(Matrix &rOp, const EPMFloat factor, const Matrix &fullOp)
   {
      rOp = factor*(fullOp.topLeftCorner(this->mNBOp, this->mNBOp) + fullOp.topRightCorner(this->mNBOp, this->mNBCs) * this->G());
   }

   void ExtensionMatrix::addRestrictedOperator(Matrix &rOp, const EPMFloat factor, const Matrix &fullOp)
   {
      rOp += factor*(fullOp.topLeftCorner(this->mNBOp, this->mNBOp) + fullOp.topRightCorner(this->mNBOp, this->mNBCs) * this->G());
   }

   void ExtensionMatrix::buildNHMatrix(const EPMFloat factor, const Matrix &fullOp)
   {
      this->mNHMatrix = factor*(fullOp.topRightCorner(this->mNBOp, this->mNBCs) * this->H());
   }

}
