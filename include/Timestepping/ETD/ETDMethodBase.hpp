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
#include "Timestepping/ETD/ETDOperator.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a ETD Metdod
    *
    * \tparam TSimType Type of the simulation
    */
   template <TSimType> class ETDMethodBase: public TimestepSchemeBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for the smart pointer to an ETDStep
         typedef EPMSHARED_PTR<ETDIteration>   SmartETDIteration;

         /// Typedef for the smart pointer to an ETDOperator
         typedef EPMSHARED_PTR<ETDOperator>   SmartETDOperator;

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
         
      protected:
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
          * @brief Do we required intermediate computations ?
          */
         bool hasIntermediate() const;

         /**
          * @brief Compute value at new time
          *
          * @param rVar Input/Output variable
          * @param nTerms New non linear terms
          */
         void doStep(ScalarType& rVar, ScalarType& nTerms);

      private:
   };

   template <typename TSimType> inline bool ETDMethodBase<TSimType>::hasIntermediate() const
   {
      return this->mCounter;
   }

   template <typename TSimType> ETDMethodBase<TSimType>::ETDMethodBase(EPMFloat a, EPMFloat b, const typename ETDMethodBase<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : TimestepSchemeBase<TSimType>(tsteps, pTrunc), mA(a), mB(b), mrBasis(basis)
   {
   }

   template <typename TSimType> void ETDMethodBase<TSimType>::doStep(typename ETDMethodBase<TSimType>::ScalarType& rVar, typename ETDMethodBase<TSimType>::ScalarType& nTerms)
   {
      this->mETDSteps.at(this->mCounter)->compute(rVar, nTerms);
   }

}

#endif // ETDMETHODBASE_HPP
