/** \file ETDMethodBase.hpp
 */

#ifndef ETDMETHODBASE_HPP
#define ETDMETHODBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/SpectralTruncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/TimestepSchemeBase.hpp"
#include "Timestepping/ETD/ETDIteration.hpp"
#include "Timestepping/ETD/ETDOperator.hpp"
#include "Simulations/SimulationTraits.hpp"

namespace EPMDynamo {

   template <TSim> class ETDMethodBase: public TimestepSchemeBase<TSim>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSim::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSim::ScalarType    ScalarType;

         /// Typedef for the smart pointer to an ETDStep
         typedef EPMDYNAMO_SHAREDPTRNS::shared_ptr<ETDIteration>   SmartETDIteration;

         /// Typedef for the smart pointer to an ETDOperator
         typedef EPMDYNAMO_SHAREDPTRNS::shared_ptr<ETDOperator>   SmartETDOperator;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          */
         ETDMethodBase(DynamoFloat a, DynamoFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartSTrunc pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDMethodBase();
         
      protected:
         /**
          * @brief \f$a\f$ coefficient of timestep operator
          */
         DynamoFloat    mA;

         /**
          * @brief \f$b\f$ coefficient of timestep operator
          */
         DynamoFloat    mB;

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
          */
         void doStep(ScalarType& rVar, ScalarType& nTerms);

      private:
   };

   template <typename TSim> inline bool ETDMethodBase<TSim>::hasIntermediate() const
   {
      return this->mCounter;
   }

   template <typename TSim> ETDMethodBase<TSim>::ETDMethodBase(DynamoFloat a, DynamoFloat b, const typename ETDMethodBase<TSim>::BasisType &basis, TimestepParameters &tsteps, SmartSTrunc pTrunc)
      : TimestepSchemeBase<TSim>(tsteps, pTrunc), mA(a), mB(b), mrBasis(basis)
   {
   }

   template <typename TSim> ETDMethodBase<TSim>::~ETDMethodBase()
   {
   }

   template <typename TSim> void ETDMethodBase<TSim>::doStep(typename ETDMethodBase<TSim>::ScalarType& rVar, typename ETDMethodBase<TSim>::ScalarType& nTerms)
   {
      this->mETDSteps.at(this->mCounter)->compute(rVar, nTerms);
   }

}

#endif // ETDMETHODBASE_HPP
