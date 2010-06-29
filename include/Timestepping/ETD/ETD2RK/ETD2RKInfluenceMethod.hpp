/** \file ETD2RKInfluenceMethod.hpp
 *  \brief Implementation of the ETD2RK method (with influence matrix)
 */

#ifndef ETD2RKINFLUENCEMETHOD_HPP
#define ETD2RKINFLUENCEMETHOD_HPP

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
#include "Timestepping/ETD/ETD2RK/ETD2RKMethod.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD2RKInfluenceMethod: public ETD2RKMethod<TSimType>
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
         ETD2RKInfluenceMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKInfluenceMethod() {};

         /**
          * @brief Add boundary condition for inlfuence matrix
          *
          * @param pBC Boundary condition
          */
         void addInfluenceBC(SmartBC pBC);

         /**
          * @brief Initialise the ETD2RK influence matrix method
          */
         void init();
         
      protected:

         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

         /**
          * @brief Compute next scheme iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         void doIteration(ScalarType& rVar, ScalarType& rNTerms);

      private:
         /**
          * @brief Number of boundary conditions required for influence matrix
          */
         int mInfluenceNBC;

         /**
          * @brief Influence matrix
          */
         InfluenceMatrix<TSimType> mInfluence;

         /**
          * @brief Initialise the influenc matrix
          */
         void initInfluence();

         /**
          * @brief update the influence matrix solution
          *
          *
          * \epmBug Not implemented yet
          */
         void updateInfluenceSolution();

         /**
          * @brief Correct timestep solution with Influence matrix solution
          */
         void useInfluenceSolution(ScalarType& rVar);
   };

   template <typename TSimType> ETD2RKInfluenceMethod<TSimType>::ETD2RKInfluenceMethod(EPMFloat a, EPMFloat b, const typename ETD2RKInfluenceMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : ETD2RKMethod<TSimType>(a, b, basis, tsteps, pTrunc), mInfluenceNBC(-1), mInfluence(pTrunc, basis)
   {
   }

   template <typename TSimType> void ETD2RKInfluenceMethod<TSimType>::addInfluenceBC(SmartBC pBC)
   {
      // Add boundary condition
      this->mInfluence.addBC(pBC);

      // Increment number of implemented influence matrix BCs
      ++this->mInfluenceNBC;
   }

   template <typename TSimType> void ETD2RKInfluenceMethod<TSimType>::init()
   {
      // Initialise the timestep operators
      ETD2RKMethod<TSimType>::init();

      // Initialise the influence matrix
      this->initInfluence();
   }

   template <typename TSimType> void ETD2RKInfluenceMethod<TSimType>::updateTimeMatrices()
   {
      // Update the timestep operators
      ETD2RKMethod<TSimType>::updateTimeMatrices();

      // Update the influence matrix solution
      this->updateInfluenceSolution();
   }

   template <typename TSimType> void ETD2RKInfluenceMethod<TSimType>::doIteration(typename ETD2RKInfluenceMethod<TSimType>::ScalarType& rVar, typename ETD2RKInfluenceMethod<TSimType>::ScalarType& rNTerms)
   {
      // Solve influence matrix part
      this->mInfluence.solve(rNTerms);

      // Go on with normal timestep
      ETDMethodBase<TSimType>::doIteration(rVar, rNTerms);

      // Correct solution with influence solution
      this->useInfluenceSolution(rVar);
   }

   template <typename TSimType> void ETD2RKInfluenceMethod<TSimType>::initInfluence()
   {
      if(this->mInfluenceNBC  == 0)
      {
         // Initialise the operators
         this->mInfluence.initOperators();

         // Compute the operators
         this->mInfluence.computeOperators();
      } else
      {
         throw EPMException("ETD2RKInfluenceMethod::initInfluence", "Tried to initialise with wrong number of BCs");
      }
   }

   template <typename TSimType> void ETD2RKInfluenceMethod<TSimType>::updateInfluenceSolution()
   {
   }

   template <typename TSimType> void ETD2RKInfluenceMethod<TSimType>::useInfluenceSolution(typename ETD2RKInfluenceMethod<TSimType>::ScalarType& rVar)
   {
      this->mInfluence.correctSolution(rVar);
   }

}

#endif // ETD2RKINFLUENCEMETHOD_HPP
