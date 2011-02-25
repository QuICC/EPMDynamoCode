/** \file ExtensionMatrix.hpp
 *  \brief Implementation of the extension matrix to impose boundary conditions
 */

#ifndef EXTENSIONMATRIX_HPP
#define EXTENSIONMATRIX_HPP

// System includes
//

// External includes
//
#include "LAPACK_Iface.hpp"

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the extension matrix to impose boundary conditions
    */
    class ExtensionMatrix
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
          * @brief Compute the missing coefficients including NH correction
          *
          * @param rVector Vector of coefficients to extend
          */
         void extend(Array& rVector) const;

         /**
          * @brief Compute the missing coefficients including NH correction
          *
          * @param rMat Matrix of coefficients to extend
          */
         void extend(MatrixZ& rMat) const;

         /**
          * @brief Compute the missing coefficients with zero BC value
          *
          * @param rVector Vector of coefficients to extend
          */
         void extendZero(Array& rVector) const;

         /**
          * @brief Compute the missing coefficients with zero BC value
          *
          * @param rMat Matrix of coefficients to extend
          */
         void extendZero(MatrixZ& rMat) const;

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

         /**
          * @brief Build the NH value imposing matrix
          *
          * @param factor Multiplicative factor
          * @param fullOp Complete operator
          */
         void buildNHMatrix(const EPMFloat factor, const Matrix& fullOp);
         
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

   inline const Matrix& ExtensionMatrix::G() const
   {
      return this->mGMatrix;
   }

   inline const Matrix& ExtensionMatrix::H() const
   {
      return this->mHMatrix;
   }

   inline const Matrix& ExtensionMatrix::NH() const
   {
      return this->mNHMatrix;
   }

}

#endif // EXTENSIONMATRIX_HPP
