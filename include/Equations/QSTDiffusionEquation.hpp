/** \file QSTDiffusionEquation.hpp
 *  \brief General QST decomposition diffusion equation
 */

#ifndef QSTDIFFUSIONEQUATION_HPP
#define QSTDIFFUSIONEQUATION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/TimeEquation.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief This class specialises the TimeEquation class to the case where the unknown
    *        is a general field (using QST decomposition)
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    */
   template <typename TSimType, typename TFieldType> class QSTDiffusionEquation : public TimeEquation<TSimType, TFieldType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local transform type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from the simulation trait to local timestepping scheme
         typedef typename SimulationTraits<TSimType>::TimestepTraits  TimeTraits;

         /**
          * @brief Constructs the Time equation from the given general Physical Field, the transform reference and the timestepping parameters
          *
          * @param rF Reference to the unknown field
          * @param transform Reference to the transform object
          * @param tsteps Timestep parameters
          * @param nBCQ Number of Q boundary conditions
          * @param nBCS Number of S boundary conditions
          * @param nBCT Number of T boundary conditions
          * @param a \f$a\f$ parameter for the timestepper
          * @param b \f$b\f$ parameter for the timestepper
          */
         QSTDiffusionEquation(TFieldType &rF, TransformType &transform, TimestepParameters &tsteps, int nBCQ, int nBCS, int nBCT, DynamoFloat a, DynamoFloat b);

         /**
          * @brief Simple empty destructor
          */
         virtual ~QSTDiffusionEquation() {};

         /**
          * @brief Timestep Equation (chooses between Predictor and corrector)
          */
         void timestep();

         /**
          * @brief Pure virtual method for initialising equation
          *
          * The implementation will be left open until a specialised equation is defined
          */
         virtual void init() = 0;

         /**
          * @brief Add a boundary condition to the Q component equation
          *
          * \param pBC Boundary condition to add to Q component list
          */
         void addQBC(SmartBC pBC);

         /**
          * @brief Add a boundary condition to the S component equation
          *
          * \param pBC Boundary condition to add to S component list
          */
         void addSBC(SmartBC pBC);

         /**
          * @brief Add a boundary condition to the T component equation
          *
          * \param pBC Boundary condition to add to T component list
          */
         void addTBC(SmartBC pBC);
         
      protected:

         /**
          * @brief Timestepper object responsible for the Q component
          */
         typename TimeTraits::Timestepper  mQTStepper;

         /**
          * @brief Timestepper object responsible for the S component
          */
         typename TimeTraits::Timestepper  mSTStepper;

         /**
          * @brief Timestepper object responsible for the T component
          */
         typename TimeTraits::Timestepper  mTTStepper;

         /**
          * @brief Vector of BoundaryCondition for the Q component
          */
         std::vector<SmartBC>    mQBCs;

         /**
          * @brief Vector of BoundaryCondition for the S component
          */
         std::vector<SmartBC>    mSBCs;

         /**
          * @brief Vector of BoundaryCondition for the T component
          */
         std::vector<SmartBC>    mTBCs;

         /**
          * @brief Initialise the timestepper objects
          */
         virtual void initTSteppers();

         /**
          * @brief  Compute the spectral expansion of the non linear terms
          */
         void transformNTerms();

      private:
   };

   template <typename TSimType, typename TFieldType> QSTDiffusionEquation<TSimType, TFieldType>::QSTDiffusionEquation(TFieldType &rF, typename QSTDiffusionEquation<TSimType, TFieldType>::TransformType &transform, TimestepParameters& tsteps, int nBCQ, int nBCS, int nBCT, DynamoFloat a, DynamoFloat b)
      : TimeEquation<TSimType, TFieldType>(rF, transform, tsteps), mQTStepper(a, b, transform.radBasis(), tsteps, rF.trunc()), mSTStepper(a, b, transform.radBasis(), tsteps, rF.trunc()), mTTStepper(a, b, transform.radBasis(), tsteps, rF.trunc())
   {
      // Set counter to total number of boundary conditions
      this->mBCCounter = nBCQ + nBCS + nBCT;
   }

   template <typename TSimType, typename TFieldType> void QSTDiffusionEquation<TSimType, TFieldType>::addQBC(SmartBC pBC)
   {
      // Add Q component BC to list
      this->mQBCs.push_back(pBC);

      // Decrement BC counter
      --this->mBCCounter;
   }

   template <typename TSimType, typename TFieldType> void QSTDiffusionEquation<TSimType, TFieldType>::addSBC(SmartBC pBC)
   {
      // Add S component BC to list
      this->mSBCs.push_back(pBC);

      // Decrement BC counter
      --this->mBCCounter;
   }

   template <typename TSimType, typename TFieldType> void QSTDiffusionEquation<TSimType, TFieldType>::addTBC(SmartBC pBC)
   {
      // Add T component BC to list
      this->mTBCs.push_back(pBC);

      // Decrement BC counter
      --this->mBCCounter;
   }

   template <typename TSimType, typename TFieldType> inline void QSTDiffusionEquation<TSimType, TFieldType>::initTSteppers()
   {
      // Set Boundary condition of timestepper for Q component
      for(int i=0; i < this->mQBCs.size(); ++i)
      {
         this->mQTStepper.addBC(this->mQBCs.at(i));
      }

      // Initialise Q component timestepper
      this->mQTStepper.init();

      // Clear list of Q component boundary conditions
      this->mQBCs.clear();


      // Set Boundary condition of timestepper for S component
      for(int i=0; i < this->mSBCs.size(); ++i)
      {
         this->mSTStepper.addBC(this->mSBCs.at(i));
      }

      // Initialise Q component timestepper
      this->mSTStepper.init();

      // Clear list of Q component boundary conditions
      this->mSBCs.clear();


      // Set Boundary condition of timestepper for T component
      for(int i=0; i < this->mTBCs.size(); ++i)
      {
         this->mTTStepper.addBC(this->mTBCs.at(i));
      }

      // Initialise Q component timestepper
      this->mTTStepper.init();

      // Clear list of Q component boundary conditions
      this->mTBCs.clear();
   }

   template <typename TSimType, typename TFieldType> inline void QSTDiffusionEquation<TSimType, TFieldType>::transformNTerms()
   {
      this->mrTransform.transformRTP2QST(this->mNTerms.rOc().rPerturbation(), this->mNTerms.oc().rtp());
   }

   template <typename TSimType, typename TFieldType> inline void QSTDiffusionEquation<TSimType, TFieldType>::timestep()
   {
      // Timestep the Q component equation
      this->mQTStepper.timestep(this->mrX.rOc().rPerturbation().rQ(), this->mNTerms.rOc().rPerturbation().rQ());

      // Timestep the S component equation
      this->mSTStepper.timestep(this->mrX.rOc().rPerturbation().rS(), this->mNTerms.rOc().rPerturbation().rS());

      // Timestep the T component equation
      this->mTTStepper.timestep(this->mrX.rOc().rPerturbation().rT(), this->mNTerms.rOc().rPerturbation().rT());
   }

}

#endif // QSTDIFFUSIONEQUATION_HPP
