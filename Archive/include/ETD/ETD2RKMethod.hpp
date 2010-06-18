/** \file ETD2RKMethod.hpp
 */

#ifndef ETD2RKMETHOD_HPP
#define ETD2RKMETHOD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/SpectralTruncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/ETDMethodBase.hpp"
#include "Timestepping/ETD/Operators/ETD2Operators.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKA.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKTimestep.hpp"
#include "Simulations/SimulationTraits.hpp"

namespace EPMDynamo {

   template <TSim> class ETD2RKMethod: public ETDMethodBase<TSim>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSim::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          */
         ETD2RKMethod(DynamoFloat a, DynamoFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartSTrunc pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKMethod();

         /**
          * @brief Reorder the previous values to use latest step
          */
         void reorderPrevious(ScalarType& rVar, ScalarType& nTerms);

         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();
         
      protected:

      private:
         /**
          * @brief Initialise the method
          */
         void initMethod();

         /**
          * @brief The set of ETD2 operators
          */
         ETD2Operators  mETD2;
   };

   template <typename TSim> ETD2RKMethod<TSim>::ETD2RKMethod(DynamoFloat a, DynamoFloat b, const typename ETD2RKMethod<TSim>::BasisType &basis, TimestepParameters &tsteps, SmartSTrunc pTrunc)
      : ETDMethodBase<TSim>(a, b, basis, tsteps, pTrunc), mETD2(pTrunc)
   {
      // Init the method
      this->initMethod();
   }

   template <typename TSim> ETD2RKMethod<TSim>::~ETD2RKMethod()
   {
   }

   template <typename TSim> void ETD2RKMethod<TSim>::initMethod()
   {
      // Create M0 operator
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<>() pOpM0;
      // Create M1 operator
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<>() pOpM1;
      // Create M2 operator
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<>() pOpM2;

      // Add required operators
         // Add M0 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(0));
         // Add M1 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(1));
         // Add M2 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(2));

      // Create storage for a variable
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<ScalarType>() pVarA;

      // Add required ETD variables
         // Add storage for variable A 
      this->mETDVars.push_back(pVarA);

      // Add required ETD NTerms
      this->mETDNTerms.push_back();

      // Create intermediate value a computation step
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<ETD2RKA>() pItA;
      // Create timestep computation step
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<ETD2RKTimestep>() pItTimestep;

      // Add required ETD steps
         // Add intermediate value A computation
      this->mETDSteps.push_back(pItA);
         // Add timestep step
      this->mETDSteps.push_back(pItTimestep);
   }

}

#endif // ETD2RKMETHOD_HPP
