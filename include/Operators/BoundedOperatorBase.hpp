/** \file BoundedOperatorBase.hpp
 *  \brief Base of the implementation of a general bounded operator
 */

#ifndef BOUNDEDOPERATORBASE_HPP
#define BOUNDEDOPERATORBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/PolyOps/POperator.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a general bounded operator (homogeneous or non homogeneous)
    *
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class BoundedOperatorBase: public TOpType
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nBC Number of boundary conditions
          * @param nTau Size of the bounded operator
          * @param id ID of the operator (for example harmonic degree l)
          */
         BoundedOperatorBase(const int nBC, const int nTau, const int id);

         /**
          * @brief Destructor
          */
         virtual ~BoundedOperatorBase() {};

         /**
          * @brief Implement the boundary conditions
          *
          * @param bcRows Matrix of the boundary imposing row values
          */
         virtual void implementBCs(const Matrix& bcRows) = 0;

         /**
          * @brief Compute the operator the bounded operator based on input
          *
          * @param factor Multiplicative factor
          * @param polyOp Polynomial operator defining the operator
          *
          * \bug PolynomialOperator needs restructuring
          */
         virtual void constructBOperator(const DynamoFloat factor, const PolynomialOperator<POperator>& polyOp) = 0; 

         /**
          * @brief Solve linear equation
          *
          * @param vector RHS of the linear equation
          */
         virtual void solve(Array& vector) = 0;
         
         /**
          * @brief Get size of the bounded operator
          */
         int nTau() const;
         
         /**
          * @brief Get identification int
          */
         int id() const;

      protected:
         /**
          * @brief Get number of boundary conditions
          */
         int nBC() const;

      private:
         /**
          * @brief Number of boundary conditions
          */
         int mNbc;

         /**
          * @brief Size of the bounded operator
          */
         int mNTau;

         /**
          * @brief Identification ID (will be and l index for example)
          */
         int mID;
   };

   template <typename TOpType> BoundedOperatorBase<TOpType>::BoundedOperatorBase(const int nBC, const int nTau, const int id)
      : TOpType(nTau), mNbc(nBC), mNTau(nTau), mID(id)
   {
   }

   template <typename TOpType> int BoundedOperatorBase<TOpType>::nBC() const
   {
      return this->mNbc;
   }

   template <typename TOpType> int BoundedOperatorBase<TOpType>::nTau() const
   {
      return this->mNTau;
   }

   template <typename TOpType> int BoundedOperatorBase<TOpType>::id() const
   {
      return this->mID;
   }

}

#endif // BOUNDEDOPERATORBASE_HPP
