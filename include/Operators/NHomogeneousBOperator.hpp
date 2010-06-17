/** \file NHomogeneousBOperator.hpp
 *  \brief Implementation of a bounded operator on a set of non homogeneous boundary conditions
 */

#ifndef NHOMOGENEOUSBOPERATOR_HPP
#define NHOMOGENEOUSBOPERATOR_HPP

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
    * \brief Implementation of a bounded operator on a set of non homogeneous boundary conditions
    *
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class NHomogeneousBOperator: public BoundedOperator<TOpType>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nBC Number of boundary conditions
          * @param size Size of the operator (basis size)
          * @param id ID of the operator (for example harmonic degree)
          */
         NHomogeneousBOperator(const int nBC, const int size, const int id);

         /**
          * @brief Destructor
          */
         virtual ~NHomogeneousBOperator() {};

         /**
          * @brief Implement the boundary conditions
          *
          * @param bcRows BC imposing rows
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

         /**
          * @brief Finalise the operator
          */
         void finaliseOp();

         /**
          * @brief Set Boundary values
          *
          * @param bcValue Value at the boundary (to impose on RHS)
          */
         virtual void setBCValue(DynamoFloat bcValue);
         
      protected:

      private:
         /**
          * @brief Boundary bordering rows values
          */
         Matrix mBCRows;

         /**
          * @brief Boundary value
          */
         DynamoFloat mBCValue;
   };

   template <typename TOpType> NHomogeneousBOperator<TOpType>::NHomogeneousBOperator(const int nBC, const int size, const int id)
      : BoundedOperator<TOpType>(nBC, size, id), mBCRows(nBC, size), mBCValue(0.0)
   {
   }

   template <typename TOpType> inline void NHomogeneousBOperator<TOpType>::implementBCs(const Matrix& bcRows)
   {
      // Simply store boundary rows
      this->mBCRows = bcRows;
   }

   template <typename TOpType> void NHomogeneousBOperator<TOpType>::constructBOperator(const DynamoFloat factor, const PolynomialOperator<POperator>& polyOp)
   {
      // Simply store polynomial operator
      this->rOp() = factor*polyOp.op();
   }

   template <typename TOpType> void NHomogeneousBOperator<TOpType>::solve(Array& vector)
   {
      // fill in value for last coefficient(s) (i.e. boundary condition)
      vector(this->nTau()-1) = this->mBCValue;

      // Call basic solve
      this->solveEquation(vector);
   }

   template <typename TOpType> void NHomogeneousBOperator<TOpType>::finaliseOp()
   {
      // Do boundary bordering
      if(this->isFactorised())
      {
         this->rOp().row(this->nTau()-this->nBC()) = this->mBCRows;
      }

      // Call the parent finalisation operation (for example LU factorisation)
      TOpType::finaliseOp();
   }

   template <typename TOpType> void NHomogeneousBOperator<TOpType>::setBCValue(DynamoFloat bcValue)
   {
      // Set non homogeneous Boundary value
      this->mBCValue = bcValue;
   }

}

#endif // NHOMOGENEOUSBOPERATOR_HPP
