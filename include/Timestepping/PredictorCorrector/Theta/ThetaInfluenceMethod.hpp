/** \file ThetaInfluenceMethod.hpp
 *  \brief Implementation of the theta method with influence matrix
 */

#ifndef THETAINFLUENCEMETHOD_HPP
#define THETAINFLUENCEMETHOD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Domain/Truncation.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaMethod.hpp"
#include "Equations/InfluenceMatrix.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method with influence matrix
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaInfluenceMethod: public ThetaMethod<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local truncation type
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
         ThetaInfluenceMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaInfluenceMethod() {};

         /**
          * @brief Add boundary condition
          *
          * The boundary conditions are added to the influence matrix. Only the first one is propagated to the standard time operators.
          *
          * @param pBC Boundary condition
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Initialise the Theta influence matrix method 
          */
         void init();
         
      protected:

         /**
          * @brief Update the timestepping operators after change of timestep length
          */
         void updateTimeMatrices();

         /**
          * @brief Set the predictor step RHS
          *
          * @param rVar Input/Ouput variable
          * @param oldVar Previous variable
          * @param rNTerms Non linear terms
          */
         void setPredictorRHS(ScalarType& rVar, const ScalarType& oldVar, ScalarType& rNTerms);

         /**
          * @brief Overload solvePredictor to include influence matrix solution step
          *
          * @param rVar Input/Ouput variable
          */
         void solvePredictor(ScalarType& rVar);

         /**
          * @brief Set the predictor step RHS
          *
          * @param rNewNTerms New non linear terms
          */
         void setCorrectorRHS(ScalarType& rNewNTerms);

         /**
          * @brief Overload solveCorrector to include influence matrix solution step
          */
         void solveCorrector();

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
          */
         void updateInfluence();
   };

   template <typename TSimType> ThetaInfluenceMethod<TSimType>::ThetaInfluenceMethod(EPMFloat a, EPMFloat b, const typename ThetaInfluenceMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : ThetaMethod<TSimType>(a, b, basis, tsteps, pTrunc), mInfluenceNBC(-2), mInfluence(pTrunc, basis)
   {
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::addBC(SmartBC pBC)
   {
      // Propagate only the first boundary condition to standard implementation
      if(this->mInfluenceNBC == -2)
      {
         ThetaMethod<TSimType>::addBC(pBC);
      }

      // Add boundary condition
      this->mInfluence.addBC(pBC);

      // Increment number of implemented influence matrix BCs
      ++this->mInfluenceNBC;
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::init()
   {
      // Initialise the timestep operators
      ThetaMethod<TSimType>::init();

      // Initialise the influence matrix
      this->initInfluence();
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::updateTimeMatrices()
   {
      // Update the timestep operators
      ThetaMethod<TSimType>::updateTimeMatrices();

      // Update the influence matrix solution
      this->updateInfluence();
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::setPredictorRHS(typename ThetaInfluenceMethod<TSimType>::ScalarType& rVar, const typename ThetaInfluenceMethod<TSimType>::ScalarType& oldVar, typename ThetaInfluenceMethod<TSimType>::ScalarType& rNTerms)
   {
      // Solve influence matrix equation
      this->mInfluence.solve(rNTerms);

      // Compute RHS part of predictor timestep equation
      ThetaMethod<TSimType>::setPredictorRHS(rVar, oldVar, rNTerms);
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::solvePredictor(typename ThetaInfluenceMethod<TSimType>::ScalarType& rVar)
   {
      // Solve for unknown variable
      ThetaMethod<TSimType>::solvePredictor(rVar);

      // Include kernel influence
      this->mInfluence.addKernel(rVar);
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::setCorrectorRHS(typename ThetaInfluenceMethod<TSimType>::ScalarType& rNewNTerms)
   {
      // Solve influenc matrix equation
      this->mInfluence.solve(rNewNTerms);
      
      // Compute RHS part of corrector timestep equation
      ThetaMethod<TSimType>::setCorrectorRHS(rNewNTerms);
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::solveCorrector()
   {
      // Get the correction to the unknown variable
      ThetaMethod<TSimType>::solveCorrector();

      // Include kernel influence
      this->mInfluence.addKernel(this->mPreviousNTerms);
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::initInfluence()
   {
      if(this->mInfluenceNBC  == 0)
      {
         // Initialise the operators
         this->mInfluence.initOperators();

         // Compute the operators
         this->mInfluence.computeOperators();
      } else
      {
         throw EPMException("ThetaInfluenceMethod::initInfluence", "Tried to initialise with wrong number of BCs");
      }
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::updateInfluence()
   {
      // Get size of radial truncation
      int nN = this->mLHS.trunc()->sim()->rad()->nN();
      // Get number of harmonic degrees
      int nL = this->mLHS.trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = ! this->mLHS.trunc()->local()->spec()->lArray()(0);

      // Create temporary storage
      Array tmp(nN);

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {
         // Initialise influence matrix solution to kernel r^l
         tmp.setConstant(0.0);
         tmp(0) = 1.0;

         // Compute the kernel influence
         this->mLHS.solveZeroVector(tmp, l);

         // Store solution from influence matrix
         this->mInfluence.storeKernelBC(tmp, l);
      }
   }

}

#endif // THETAINFLUENCEMETHOD_HPP
