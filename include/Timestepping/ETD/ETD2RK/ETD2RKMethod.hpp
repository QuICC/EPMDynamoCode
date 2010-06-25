/** \file ETD2RKMethod.hpp
 *  \brief Implementation of the ETD2RK method (without influence matrix)
 */

#ifndef ETD2RKMETHOD_HPP
#define ETD2RKMETHOD_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/ETDMethodBase.hpp"
#include "Timestepping/ETD/Operators/ETD2Operators.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKA.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKTimestep.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    *
    * \tparam TSimType Type of the simulation
    */
   template <TSimType> class ETD2RKMethod: public ETDMethodBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          */
         ETD2RKMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKMethod() {};

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

   template <typename TSimType> ETD2RKMethod<TSimType>::ETD2RKMethod(EPMFloat a, EPMFloat b, const typename ETD2RKMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : ETDMethodBase<TSimType>(a, b, basis, tsteps, pTrunc), mETD2(pTrunc)
   {
      // Init the method
      this->initMethod();
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::initMethod()
   {
      // Create M0 operator
      EPMSHARED_PTR<>() pOpM0;
      // Create M1 operator
      EPMSHARED_PTR<>() pOpM1;
      // Create M2 operator
      EPMSHARED_PTR<>() pOpM2;

      // Add required operators
         // Add M0 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(0));
         // Add M1 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(1));
         // Add M2 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(2));

      // Create storage for a variable
      EPMSHARED_PTR<ScalarType>() pVarA;

      // Add required ETD variables
         // Add storage for variable A 
      this->mETDVars.push_back(pVarA);

      // Add required ETD NTerms
      this->mETDNTerms.push_back();

      // Create intermediate value a computation step
      EPMSHARED_PTR<ETD2RKA>() pItA;
      // Create timestep computation step
      EPMSHARED_PTR<ETD2RKTimestep>() pItTimestep;

      // Add required ETD steps
         // Add intermediate value A computation
      this->mETDSteps.push_back(pItA);
         // Add timestep step
      this->mETDSteps.push_back(pItTimestep);
   }

}

#endif // ETD2RKMETHOD_HPP
