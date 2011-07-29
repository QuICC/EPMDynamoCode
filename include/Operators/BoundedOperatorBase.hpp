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
#include "BoundaryConditions/BoundaryCondition.hpp"

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
          * @param simId Simulation wide ID of the operator (for example harmonic degree l)
          */
         BoundedOperatorBase(const int nBC, const int nTau, const int id, const int simID);

         /**
          * @brief Destructor
          */
         virtual ~BoundedOperatorBase() {};

         /**
          * @brief Implement the boundary conditions into operator
          */
         virtual void implementBCs() = 0;

         /**
          * @brief Add a boundary condition to operator
          *
          * @param pBC Smart pointer to boundary condition
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Set operator's diagonal value
          */
         virtual void setOperator(const EPMFloat factor) = 0;

         /**
          * @brief Set operator value
          */
         virtual void setOperator(const EPMFloat factor, const Matrix& mat) = 0;

         /**
          * @brief Add operator's diagonal value
          */
         virtual void addOperator(const EPMFloat factor) = 0;

         /**
          * @brief Add operator value
          */
         virtual void addOperator(const EPMFloat factor, const Matrix& mat) = 0;

         /**
          * @brief Restrict the operator by including the boundary conditions
          */
         virtual void constructBoundedOperator() = 0; 

         /**
          * @brief Solve linear equation for a homogeneous BC
          *
          * @param vector RHS of the linear equation
          */
         virtual void solveZero(Array& vector) = 0;

         /**
          * @brief Solve linear equation
          *
          * @param vector RHS of the linear equation
          * @param isReal Is real component of equation?
          */
         virtual void solve(Array& vector, const int simM, const bool isReal) = 0;
         
         /**
          * @brief Get size of the bounded operator
          */
         int nTau() const;
         
         /**
          * @brief Get identification int
          */
         int id() const;
         
         /**
          * @brief Simulation wide identification id
          */
         int simId() const;

      protected:
         /**
          * @brief Get number of boundary conditions
          */
         int nBC() const;

         /**
          * @brief Vector of BoundaryCondition
          */
         std::vector<SmartBC>    mBCs;

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

         /**
          * @brief Simulation wide identification ID (will be and l index for example)
          */
         int mSimID;
   };

   template <typename TOpType> BoundedOperatorBase<TOpType>::BoundedOperatorBase(const int nBC, const int nTau, const int id, const int simID)
      : TOpType(nTau), mNbc(nBC), mNTau(nTau), mID(id), mSimID(simID)
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

   template <typename TOpType> int BoundedOperatorBase<TOpType>::simId() const
   {
      return this->mSimID;
   }

   template <typename TOpType> void BoundedOperatorBase<TOpType>::addBC(SmartBC pBC) 
   {
      // Add BC
      this->mBCs.push_back(pBC);
   }

}

#endif // BOUNDEDOPERATORBASE_HPP
