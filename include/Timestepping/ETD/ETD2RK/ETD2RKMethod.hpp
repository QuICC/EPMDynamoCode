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
#include "Timestepping/IterativeSchemeBase.hpp"
#include "Timestepping/ETD/Operators/ETD2Operators.hpp"
#include "Timestepping/ETD/Iterations/ETD1Iteration.hpp"
#include "Timestepping/ETD/Iterations/ETD2RKTimestep.hpp"
#include "Config/SimulationConfig.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD2RKMethod: public IterativeSchemeBase<TSimType>
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
          * @param hasL0 Is l=0 mode required?
          */
         ETD2RKMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

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
          * @brief The set of ETD2 operators
          */
         ETD2Operators<TSimType>  mETD2;

         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

         /**
          * @brief Init the operators
          */
         void initOperators();

      private:

         /**
          * @brief Initialise the storage for the method
          */
         void initStorage();
   };

   template <typename TSimType> ETD2RKMethod<TSimType>::ETD2RKMethod(EPMFloat a, EPMFloat b, const typename ETD2RKMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : IterativeSchemeBase<TSimType>(a, b, basis, tsteps, pTrunc, hasL0), mETD2(b/a, pTrunc, hasL0)
   {
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::addBC(SmartBC pBC)
   {
      // add boundary condition to operators
      this->mETD2.addBC(pBC);
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::init()
   {
      // initialise pointers
      this->initStorage();

      // Initialise the operators
      this->initOperators();
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::initOperators()
   {
      // initialise the operators
      this->mETD2.initOperators(this->mrBasis);
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::updateTimeMatrices()
   {
      // Update the time depended matrices
      this->mETD2.update(this->rTSParams().dt(), this->mrBasis);
   }

   template <typename TSimType> void ETD2RKMethod<TSimType>::initStorage()
   {
      // Create intermediate value a computation step
      EPMSHARED_PTR<ETD1Iteration<TSimType> > pItA(new ETD1Iteration<TSimType> (1.0/this->mA, this->mETD2.pEtdF(0), this->mETD2.pEtdF(1)));

      // Create timestep computation step
      EPMSHARED_PTR<ETD2RKTimestep<TSimType> > pItTimestep(new ETD2RKTimestep<TSimType> (1.0/this->mA, this->pOldNTerms(), this->mETD2.pEtdF(2)));

      // Add required ETD steps
         // Add intermediate value A computation
      this->mSteps.push_back(pItA);
         // Add timestep step
      this->mSteps.push_back(pItTimestep);
   }

}

#endif // ETD2RKMETHOD_HPP
