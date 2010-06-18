/** \file DenseLUOperator.cpp
 *  \brief Source of the dense LU factorisable operator
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "Operators/DenseLUOperator.hpp"
#include "LAPACK_Iface.hpp"

// Project includes
//

namespace EPMDynamo {

   DenseLUOperator::DenseLUOperator(const int size)
      : DenseOperator(size, true), mIpiv(size)
   {
   }

   void DenseLUOperator::factorise()
   {
      int rows = this->op().rows();
      int cols = this->op().cols();

      int info;

      // Call LAPACK dgetrf routine for factorisation
      dgetrf_(&rows, &cols, this->rOp().data(), &rows, this->rIpiv().data(), &info);

      // Test success of computation through assert
      assert(info == 0);
   }

   void DenseLUOperator::inverse()
   {
      int n = this->op().rows();

      int info;
      int lwork = -1;
      Array work(n);

      // Call LAPACK dgetri routine for optimal work size
      dgetri_(&n, this->rOp().data(), &n, this->rIpiv().data(), work.data(), &lwork, &info);

      // Test success of computation through assert
      assert(info == 0);

      // Set to optimal work size
      lwork = static_cast<int>(work(0));
      work.resize(lwork);

      // Call LAPACK dgetri routine for inverse computation
      dgetri_(&n, this->rOp().data(), &n, this->rIpiv().data(), work.data(), &lwork, &info);

      // Test success of computation through assert
      assert(info == 0);
   }

   void DenseLUOperator::solveEquation(Array& vector)
   {
      // Type of system for LAPACK routine
      char trans = 'N';

      // Size of the linear problem
      int rows = this->op().rows();
      // Lapack info output
      int info;
      // Number of columns to solve for
      int nCol = 1;

      // Call LAPACK dgetrs to solve real part equation
      dgetrs_(&trans, &rows, &nCol, this->rOp().data(), &rows, this->rIpiv().data(), vector.data(), &rows, &info);

      // Test success of computation through assert
      assert(info == 0);
   }

}
