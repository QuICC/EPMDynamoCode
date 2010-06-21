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
      this->mHMatrix = bcRows.corner(Eigen::BottomLeft, this->mNBCs, this->mNBCs);

      // Fill in the values for G
      this->mGMatrix = bcRows.corner(Eigen::TopLeft, this->mNBCs, this->mNBOp);

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
      rVector.start(this->mNBOp) -= this->NH() * rVector.end(this->mNBCs);
   }

   void ExtensionMatrix::extendZero(Array &rVector)
   {
      rVector.end(this->mNBCs) = this->G() * rVector.start(this->mNBOp);
   }

   void ExtensionMatrix::extend(Array &rVector)
   {
      rVector.end(this->mNBCs) = this->G() * rVector.start(this->mNBOp) + this->H() * rVector.end(this->mNBCs);
   }

   void ExtensionMatrix::restrictOperator(Matrix &rOp, const EPMFloat factor, const Matrix &fullOp)
   {
      rOp = factor*(fullOp.corner(Eigen::TopLeft, this->mNBOp, this->mNBOp) + fullOp.corner(Eigen::TopRight, this->mNBOp, this->mNBCs) * this->G());
   }

   void ExtensionMatrix::buildNHMatrix(const EPMFloat factor, const Matrix &fullOp)
   {
      this->mNHMatrix = factor*(fullOp.corner(Eigen::TopRight, this->mNBOp, this->mNBCs) * this->H());
   }

}
