/** \file ThetaInfluenceMethod.hpp
 *  \brief Implementation of the predictor-corrector theta method with influence matrix
 */

#ifndef THETAINFLUENCEMETHOD_HPP
#define THETAINFLUENCEMETHOD_HPP

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
#include "Timestepping/PredictorCorrector/Theta/ThetaMethod.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluencePredictor.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluenceCorrector.hpp"
#include "Equations/InfluenceMatrix.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the predictor-corrector theta method with influence matrix
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaInfluenceMethod: public ThetaMethod<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for a smart pointer to an influence matrix operator
         typedef EPMSHARED_PTR<InfluenceMatrix<TSimType> > SmartInfluenceMatrix;

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
         ThetaInfluenceMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ThetaInfluenceMethod() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC);

         /**
          * @brief Initialise the method
          */
         void init();
         
      protected:
         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

      private:
         /**
          * @brief Number of boundary conditions required for influence matrix
          */
         int mInfluenceNBC;

         /**
          * @brief Influence matrix
          */
         SmartInfluenceMatrix mpInfluence;

         /**
          * @brief Initialise the storage for the method
          */
         void initStorage();

         /**
          * @brief Initialise the influenc matrix
          */
         void initInfluence();

         /**
          * @brief update the influence matrix solution
          */
         void updateInfluence();
   };

   template <typename TSimType> ThetaInfluenceMethod<TSimType>::ThetaInfluenceMethod(EPMFloat a, EPMFloat b, const typename ThetaInfluenceMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : ThetaMethod<TSimType>(a, b, basis, tsteps, pTrunc, hasL0), mInfluenceNBC(-2), mpInfluence(new InfluenceMatrix<TSimType> (pTrunc, basis, hasL0))
   {
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::addBC(SmartBC pBC)
   {
      // Propagate only the first boundary condition to standard implementation
      if(this->mInfluenceNBC == -2)
      {
         ThetaMethod<TSimType>::addBC(pBC);
      }

      // Add boundary condition to influence matrix
      this->mpInfluence->addBC(pBC);

      // Increment number of implemented influence matrix BCs
      ++this->mInfluenceNBC;
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::init()
   {
      // initialise pointers
      this->initStorage();

      // Initialise the operators
      this->initOperators();

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

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::initStorage()
   {
      //
      // Create intermediate storage
      EPMSHARED_PTR<ScalarType > pPrevious(new ScalarType (this->mpLHS->trunc(), this->mpLHS->hasL0()));
      this->mTmpScalars.push_back(pPrevious);

      // Create the predictor computation step
      EPMSHARED_PTR<ThetaInfluencePredictor<TSimType> > pItP(new ThetaInfluencePredictor<TSimType> (this->mpInfluence, this->mpLHS, this->mpRHS));

      // Create the corrector computation step
      EPMSHARED_PTR<ThetaInfluenceCorrector<TSimType> > pItC(new ThetaInfluenceCorrector<TSimType> (this->mpInfluence, pPrevious, this->mpLHS, this->mpRHS));

      //
      // Add the required computation steps
      
      // Add the predictor computation
      this->mSteps.push_back(pItP);
      // Add a corrector step
      this->mSteps.push_back(pItC);
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::initInfluence()
   {
      if(this->mInfluenceNBC  == 0)
      {
         // Initialise the operators
         this->mpInfluence->initOperators();

         // Compute the operators
         this->mpInfluence->computeOperators();
      } else
      {
         throw EPMException("ThetaInfluenceMethod::initInfluence", "Tried to initialise with wrong number of BCs");
      }
   }

   template <typename TSimType> void ThetaInfluenceMethod<TSimType>::updateInfluence()
   {
      // Get size of radial truncation
      int nN = this->mpLHS->trunc()->sim()->rad()->nN();
      // Get number of harmonic degrees
      int nL = this->mpLHS->trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = ! this->mpLHS->trunc()->local()->spec()->lArray()(0);

      // Create temporary storage
      Array tmp(nN);

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {
         // Initialise influence matrix solution to kernel r^l
         tmp.setConstant(0.0);
         tmp(0) = 1.0;

         // Compute the kernel influence
         this->mpLHS->solveZeroVector(tmp, l);

         // Store solution from influence matrix
         this->mpInfluence->storeKernelBC(tmp, l);
      }
   }

}

#endif // THETAINFLUENCEMETHOD_HPP
