/** \file ExtensionMatrix.hpp
 *  \brief Implementation of the extension matrix to impose boundary conditions
 */

#ifndef EXTENSIONMATRIX_HPP
#define EXTENSIONMATRIX_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the extension matrix to impose boundary conditions
    *
    * \tparam TBCType Type of the boundary condition
    */
   template <BCType TBCType> class ExtensionMatrix
   {
      public:
         /**
          * @brief Constructor
          */
         ExtensionMatrix();

         /**
          * @brief Destructor
          */
         virtual ~ExtensionMatrix() {};

         /**
          * @brief Build the extension matrices
          *
          * @param bcRows Matrix of boundary condition values
          */
         void build(const Matrix& bcRows);

         /**
          * @brief Compute the missing coefficients
          *
          * @param rVector Vector of coefficients to extend
          */
         void extend(Array& rVector);

         /**
          * @brief Compute the missing coefficients with zero BC value
          *
          * @param rVector Vector of coefficients to extend
          */
         void extendZero(Array& rVector);

         /**
          * @brief Impose boundary value on RHS
          *
          * @param rVector RHS vector
          */
         void prepareRHS(Array& rVector);

         /**
          * @brief Restrict operator
          *
          * @param rOp Restricted operator storage
          * @param factor Multiplicative factor
          * @param fullOp Complete operator
          */
         void restrictOperator(Matrix &rOp, const EPMFloat factor, const Matrix& fullOp);
         
      protected:

         /**
          * @brief Get the matrix G
          */
         const Matrix& G() const;

         /**
          * @brief Get the matrix H
          */
         const Matrix& H() const;

         /**
          * @brief Get the matrix NH
          */
         const Matrix& NH() const;

      private:
         /**
          * @brief Number of BCs
          */
         int   mNBCs;

         /**
          * @brief Size of bounded operator
          */
         int mNBOp;

         /**
          * @brief G matrix to recover remove coefficients
          */
         Matrix   mGMatrix;

         /**
          * @brief H matrix to handle the inhomogeneous boundary conditions
          */
         Matrix   mHMatrix;

         /**
          * @brief H matrix to handle the inhomogeneous boundary conditions
          */
         Matrix   mNHMatrix;

         /**
          * @brief Invert the H matrix
          */
         void invertHMatrix();
   };

   template <BCType TBCType> inline const Matrix& ExtensionMatrix<TBCType>::G() const
   {
      return this->mGMatrix;
   }

   template <BCType TBCType> inline const Matrix& ExtensionMatrix<TBCType>::H() const
   {
      return this->mHMatrix;
   }

   template <BCType TBCType> inline const Matrix& ExtensionMatrix<TBCType>::NH() const
   {
      return this->mNHMatrix;
   }

   template <BCType TBCType> ExtensionMatrix<TBCType>::ExtensionMatrix()
      : mNBCs(0), mNBOp(0)
   {
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::build(const Matrix &bcRows)
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

   template <BCType TBCType> void ExtensionMatrix<TBCType>::invertHMatrix()
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

   template <BCType TBCType> void ExtensionMatrix<TBCType>::prepareRHS(Array &rVector)
   {
      if(TBCType == NonHomogeneous)
      {
         rVector.start(this->mNBOp) -= this->NH() * rVector.end(this->mNBCs);
      }
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::extendZero(Array &rVector)
   {
      rVector.end(this->mNBCs) = this->G() * rVector.start(this->mNBOp);
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::extend(Array &rVector)
   {
      if(TBCType == Homogeneous)
      {
         rVector.end(this->mNBCs) = this->G() * rVector.start(this->mNBOp);
      } else
      {
         rVector.end(this->mNBCs) = this->G() * rVector.start(this->mNBOp) + this->H() * rVector.end(this->mNBCs);
      }
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::restrictOperator(Matrix &rOp, const EPMFloat factor, const Matrix &fullOp)
   {
      rOp = factor*(fullOp.corner(Eigen::TopLeft, this->mNBOp, this->mNBOp) + fullOp.corner(Eigen::TopRight, this->mNBOp, this->mNBCs) * this->G());

      if(TBCType == NonHomogeneous)
      {
         this->mNHMatrix = fullOp.corner(Eigen::TopRight, this->mNBOp, this->mNBCs) * this->H();
      }
   }

}

#endif // EXTENSIONMATRIX_HPP
