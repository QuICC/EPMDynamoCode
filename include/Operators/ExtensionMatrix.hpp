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
    *
    * \bug Computations are not yet implemented
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

      private:
         /**
          * @brief G matrix to recover remove coefficients
          */
         Matrix   mGMatrix;

         /**
          * @brief H matrix to handle the inhomogeneous boundary conditions
          */
         Matrix mHMatrix;
   };

   template <BCType TBCType> ExtensionMatrix<TBCType>::ExtensionMatrix()
   {
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::build(const Matrix &bcRows)
   {
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::prepareRHS(Array &rVector)
   {
      if(TBCType == NonHomogeneous)
      {
      }
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::extend(Array &rVector)
   {
   }

   template <BCType TBCType> void ExtensionMatrix<TBCType>::restrictOperator(Matrix &rOp, const EPMFloat factor, const Matrix &fullOp)
   {
   }

}

#endif // EXTENSIONMATRIX_HPP
