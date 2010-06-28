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
   template <typename TSimType> class ETD2RKMethod: public ETDMethodBase<TSimType>
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
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC);

         /**
          * @brief Initialise the ETD2RK method
          */
         void init();
         
      protected:

         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

      private:

         /**
          * @brief The set of ETD2 operators
          */
         ETD2Operators<TSimType>  mETD2;

         /**
          * @brief Initialise the storage for the method
          */
         void initStorage(SmartTruncation pTrunc);
   };

   template <typename TSimType> ETD2RKMethod<TSimType>::ETD2RKMethod(EPMFloat a, EPMFloat b, const typename ETD2RKMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : ETDMethodBase<TSimType>(a, b, basis, tsteps, pTrunc), mETD2(pTrunc)
   {
      // initialise pointers
      this->initStorage(pTrunc);
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::addBC(SmartBC pBC)
   {
      // add boundary condition to operators
      this->mETD2.addBC(pBC);
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::init()
   {
      // initialise the operators
      this->mETD2.initOperators();
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::updateTimeMatrices()
   {
      // Update the time depended matrices
      this->mETD2.update(this->rTSParams().dt(), this->mrBasis);
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::initStorage(SmartTruncation pTrunc)
   {
      // Add required operators
         // Add M0 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(0));
         // Add M1 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(1));
         // Add M2 operator
      this->mETDOperators.push_back(this->mETD2.pEtdF(2));

      // Create storage for a variable
      EPMSHARED_PTR<ScalarType> pVarA(new ScalarType(pTrunc));

      // Add required ETD variables
         // Add storage for variable A 
      this->mETDVars.push_back(pVarA);

      // Create intermediate value a computation step
      EPMSHARED_PTR<ETD2RKA<TSimType> > pItA(new ETD2RKA<TSimType> (this->mETDOperators.at(0), this->mETDOperators.at(1)));
      // Create timestep computation step
      EPMSHARED_PTR<ETD2RKTimestep<TSimType> > pItTimestep(new ETD2RKTimestep<TSimType> (this->mETDVars.at(0), this->mETDOperators.at(2)));

      // Add required ETD steps
         // Add intermediate value A computation
      this->mETDSteps.push_back(pItA);
         // Add timestep step
      this->mETDSteps.push_back(pItTimestep);
   }

}

#endif // ETD2RKMETHOD_HPP
