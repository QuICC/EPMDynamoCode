/** \file ETD1Method.hpp
 *  \brief Implementation of the ETD1 method (without influence matrix)
 */

#ifndef ETD1METHOD_HPP
#define ETD1METHOD_HPP

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
#include "Timestepping/ETD/Operators/ETD1Operators.hpp"
#include "Timestepping/ETD/Iterations/ETD1Iteration.hpp"
#include "Config/SimulationConfig.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD1Method: public IterativeSchemeBase<TSimType>
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
         ETD1Method(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETD1Method() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC);

         /**
          * @brief Initialise the ETD1 method
          */
         void init();
         
      protected:

         /**
          * @brief The set of ETD1 operators
          */
         ETD1Operators<TSimType>  mETD1;

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

   template <typename TSimType> ETD1Method<TSimType>::ETD1Method(EPMFloat a, EPMFloat b, const typename ETD1Method<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : IterativeSchemeBase<TSimType>(a, b, basis, tsteps, pTrunc, hasL0), mETD1(b/a, pTrunc, hasL0)
   {
   }

   template <typename TSimType> void ETD1Method<TSimType>::addBC(SmartBC pBC)
   {
      // add boundary condition to operators
      this->mETD1.addBC(pBC);
   }

   template <typename TSimType> void ETD1Method<TSimType>::init()
   {
      // initialise pointers
      this->initStorage();

      // Initialise the operators
      this->initOperators();
   }

   template <typename TSimType> void ETD1Method<TSimType>::initOperators()
   {
      // initialise the operators
      this->mETD1.initOperators(this->mrBasis);
   }

   template <typename TSimType> void ETD1Method<TSimType>::updateTimeMatrices()
   {
      // Update the time depended matrices
      this->mETD1.update(this->rTSParams().dt(), this->mrBasis);
   }

   template <typename TSimType> void ETD1Method<TSimType>::initStorage()
   {
      // Create intermediate value a computation step
      EPMSHARED_PTR<ETD1Iteration<TSimType> > pIt(new ETD1Iteration<TSimType> (1.0/this->mA, this->mETD1.pEtdF(0), this->mETD1.pEtdF(1)));

      // Add required ETD steps
         // Add intermediate value A computation
      this->mSteps.push_back(pIt);
   }

}

#endif // ETD1METHOD_HPP
