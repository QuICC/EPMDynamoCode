/** \file ETDMethodBase.hpp
 *  \brief Base of the implementation of a ETD Metdod
 */

#ifndef ETDMETHODBASE_HPP
#define ETDMETHODBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/TimestepSchemeBase.hpp"
#include "Timestepping/ETD/ETDIteration.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a ETD Metdod
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETDMethodBase: public TimestepSchemeBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for a smart pointer to a scalar type
         typedef EPMSHARED_PTR<ScalarType>   SmartScalarType;

         /// Typedef for the smart pointer to an ETDStep
         typedef EPMSHARED_PTR<ETDIteration<TSimType> >   SmartETDIteration;

         /// Typedef for the smart pointer to an ETDOperator
         typedef EPMSHARED_PTR<typename ETDSchemeTraits<TSimType>::Operators>   SmartETDOperator;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          */
         ETDMethodBase(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDMethodBase() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC) = 0;
         
      protected:
         /**
          * @brief  Iterator to current iteration
          */
         typename std::vector<SmartETDIteration>::iterator  mCurrentIt;

         /**
          * @brief \f$a\f$ coefficient of timestep operator
          */
         EPMFloat    mA;

         /**
          * @brief \f$b\f$ coefficient of timestep operator
          */
         EPMFloat    mB;

         /**
          * @brief Reference to the radial basis
          */
         const BasisType&   mrBasis;

         /**
          * @brief Vector of all the steps required for the selected scheme
          */
         std::vector<SmartETDIteration>  mETDSteps;

         /**
          * @brief Vector of operators required by the ETD scheme
          */
         std::vector<SmartETDOperator>  mETDOperators;

         /**
          * @brief Storage for the required variable values
          */
         std::vector<SmartScalarType>  mETDVars;

         /**
          * @brief Storage for the required NTerms
          */
         std::vector<SmartScalarType>  mETDNTerms;

         /**
          * @brief Reset current pointer to first step
          */
         void resetIterations();

         /**
          * @brief Compute next scheme iteration
          *
          * @param rVar Input/Output variable
          * @param nTerms New non linear terms
          */
         void doIteration(ScalarType& rVar, ScalarType& nTerms);

      private:
   };

   template <typename TSimType> ETDMethodBase<TSimType>::ETDMethodBase(EPMFloat a, EPMFloat b, const typename ETDMethodBase<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : TimestepSchemeBase<TSimType>(tsteps, pTrunc), mA(a), mB(b), mrBasis(basis)
   {
   }

   template <typename TSimType> void ETDMethodBase<TSimType>::resetIterations()
   {
      this->mCurrentIt = this->mETDSteps.begin();
   }

   template <typename TSimType> void ETDMethodBase<TSimType>::doIteration(typename ETDMethodBase<TSimType>::ScalarType& rVar, typename ETDMethodBase<TSimType>::ScalarType& nTerms)
   {
      // Do step computation
      (*this->mCurrentIt)->compute(rVar, nTerms);

      // Go forward one step
      ++this->mCurrentIt;
   }

}

#endif // ETDMETHODBASE_HPP
