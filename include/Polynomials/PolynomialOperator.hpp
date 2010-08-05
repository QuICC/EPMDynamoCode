/** \file PolynomialOperator.hpp
 *  \brief Implementation of a general polynomial operator
 */

#ifndef POLYNOMIALOPERATOR_HPP
#define POLYNOMIALOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/PolyOps/OperatorStorageTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a general polynomial operator
    *
    * \tparam TExp Type of the expression
    *
    * \epmTodo This NEEDS a big thought about a better and cleaner way to implement the required features!
    */
   template <typename TExp> class PolynomialOperator : public OperatorStorageTraits<TExp>::StorageType
   {
      public:
         /**
         * @brief Constructor for a P operator
         *
         * @param poly Matrix operator
         */
         PolynomialOperator(const Matrix &poly);

         /**
         * @brief Constructor for a C P operator
         *
         * @param c Multiplicative constant
         * @param poly Matrix operator
         */
         PolynomialOperator(const typename TExp::ConstantType c, const Matrix &poly);

         /**
         * @brief Constructor for a RP or LP operator
         *
         * @param arr LHS or RHS diagonal multiplication
         * @param poly Matrix operator
         */
         PolynomialOperator(const Array &arr, const Matrix &poly);

         /**
         * @brief Constructor for a CRP or CLP operator
         *
         * @param c Multiplicative constant
         * @param arr LHS or RHS diagonal multiplication
         * @param poly Matrix operator
         */
         PolynomialOperator(const typename TExp::ConstantType c, const Array &arr, const Matrix &poly);

         /**
         * @brief Constructor for a LRP operator
         *
         * @param left LHS diagonal multiplication
         * @param right RHS diagonal multiplication
         * @param poly Matrix operator
         */
         PolynomialOperator(const Array &left, const Array &right, const Matrix &poly);

         /**
         * @brief Constructor for a CLRP operator
         *
         * @param c Multiplicative constant
         * @param left LHS diagonal multiplication
         * @param right RHS diagonal multiplication
         * @param poly Matrix operator
         */
         PolynomialOperator(const typename TExp::ConstantType c, const Array &left, const Array &right, const Matrix &poly);

         /**
         * @brief Destructor
         */
         virtual ~PolynomialOperator() {};
         
      protected:
         /**
         * @brief Constructor
         */
         PolynomialOperator();

      private:

   };

   template <typename TExp> PolynomialOperator<TExp>::PolynomialOperator(const Matrix &poly)
      : OperatorStorageTraits<TExp>::StorageType(poly)
   {
   }

   template <typename TExp> PolynomialOperator<TExp>::PolynomialOperator(const typename TExp::ConstantType c, const Matrix &poly)
      : OperatorStorageTraits<TExp>::StorageType(c, poly)
   {
   }

   template <typename TExp> PolynomialOperator<TExp>::PolynomialOperator(const Array &arr, const Matrix &poly)
      : OperatorStorageTraits<TExp>::StorageType(arr, poly)
   {
   }

   template <typename TExp> PolynomialOperator<TExp>::PolynomialOperator(const typename TExp::ConstantType c, const Array &arr, const Matrix &poly)
      : OperatorStorageTraits<TExp>::StorageType(c, arr, poly)
   {
   }

   template <typename TExp> PolynomialOperator<TExp>::PolynomialOperator(const Array &left, const Array &right, const Matrix &poly)
      : OperatorStorageTraits<TExp>::StorageType(left, right, poly)
   {
   }

   template <typename TExp> PolynomialOperator<TExp>::PolynomialOperator(const typename TExp::ConstantType c, const Array &left, const Array &right, const Matrix &poly)
      : OperatorStorageTraits<TExp>::StorageType(c, left, right, poly)
   {
   }

}

#endif // POLYNOMIALOPERATOR_HPP
