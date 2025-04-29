/** \file QuadratureRule.cpp
 *  \brief Source of the implementation of a general quadrature rule computation
 */

// System includes
//

// External includes
//

// Interfaces includes
//
#include "LAPACK_Iface.hpp"

// Class include
//
#include "Polynomials/QuadratureRule.hpp"

// Project includes
//

namespace EPMDynamo {

   QuadratureRule::QuadratureRule(const Array& diagonal, const Array& subdiagonal, EPMFloat intWeight)
      :mDiagonal(diagonal), mSubdiagonal(subdiagonal), mIntWeight(intWeight)
   {
   }

   void QuadratureRule::computeGrid(Array& grid)
   {
      // LAPACK dstevx inputs
      //
      // Compute all eigenvalues and eigenvectors
      char  jobz = 'N';
      char  range = 'A';
      char  prec = 'S';

      // Sizes of the matrix
      int   n = mDiagonal.size();
      int   ldz = mDiagonal.size();
      // Get maximal precision
      EPMFloat   abstol = 2.0*dlamch_(&prec);

      // info
      int            info, m;

      // unused variables (used if given range for eigenvalues/vectors)
      int            ilu;
      EPMFloat    vlu;


      // Required arrays and matrices
      ArrayI   iFail(n);
      ArrayI   iWork(5*n);
      Array    w(n);
      Array    work(5*n);
      Array    z;

      // Call to LAPACK double real symmetric eigenvalue solver DSTEVX
      dstevx_(&jobz, &range, &n, mDiagonal.data(), mSubdiagonal.data(), &vlu, &vlu, &ilu, &ilu, &abstol, &m, w.data(), z.data(), &ldz, work.data(), iWork.data(), iFail.data(), &info);

      // Check computation with assert
      assert(info == 0);

      if(grid.size() == w.size())
      {
         grid = w;
      } else
      {
         grid = w.head(grid.size());
      }
   }

   void QuadratureRule::computeGridAndWeights(Array& grid, Array& weights)
   {
      // LAPACK dstevx inputs
      //
      // Compute all eigenvalues and eigenvectors
      char  jobz = 'V';
      char  range = 'A';
      char  prec = 'S';

      // Sizes of the matrix
      int   n = mDiagonal.size();
      int   ldz = mDiagonal.size();
      // Get maximal precision
      EPMFloat   abstol = 2.0*dlamch_(&prec);

      // info and unused variables (used if given range for eigenvalues/vectors)
      int            info, m;

      int            ilu;
      EPMFloat    vlu;


      // Required arrays and matrices
      ArrayI   iFail(n);
      ArrayI   iWork(5*n);
      Array    w(n);
      Array    work(5*n);
      Matrix    z(n,n);

      // Call to LAPACK double real symmetric eigenvalue solver DSTEVX
      dstevx_(&jobz, &range, &n, mDiagonal.data(), mSubdiagonal.data(), &vlu, &vlu, &ilu, &ilu, &abstol, &m, w.data(), z.data(), &ldz, work.data(), iWork.data(), iFail.data(), &info);

      // Check computation with assert
      assert(info == 0);

      if(grid.size() == w.size())
      {
         grid = w;
         weights = z.row(0).array() * z.row(0).array() * this->mIntWeight;
      } else
      {
         grid = w.head(grid.size());
         weights = z.row(0).head(grid.size()).array() * z.row(0).head(grid.size()).array() * this->mIntWeight;
      }
   }

}
