/** \file TorPolDiffusionEquation.hpp
 *  \brief General Toroidal/Poloidal decomposition diffusion equation
 */

#ifndef TORPOLDIFFUSIONEQUATION_HPP
#define TORPOLDIFFUSIONEQUATION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "GeneralFields/TorPolField.hpp"
#include "Equations/TimeEquation.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/SimpleTTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief This is a specialisation of the TimeEquation class for a Toroidal/Poloidal
    *        unknown field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    * \tparam TInfluenceTraits Traits for the influence matrix step
    */
   template <typename TSimType, typename TFieldType, template <typename> class TInfluenceTraits = SimpleTTraits> class TorPolDiffusionEquation : public TimeEquation<TSimType, TFieldType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local transform type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local timestepper type
         typedef typename SimulationTraits<TSimType>::TimestepTraits    TimeTraits;

         /**
          * @brief Constructs the timestepable equation 
          *
          * @param rF Unknown Divergence free field
          * @param transform Reference to the transform object
          * \param tsteps Timestep parameters
          * @param nBCT Number of toroidal boundary conditions
          * @param nBCP Number of poloidal boundary conditions
          * @param a \f$a\f$ parameter of the Timestepper
          * @param b \f$a\f$ parameter of the Timestepper
          */
         TorPolDiffusionEquation(TFieldType &rF, TransformType &transform, TimestepParameters &tsteps, int nBCT, int nBCP, EPMFloat a, EPMFloat b);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TorPolDiffusionEquation() {};

         /**
          * @brief Timestep equation
          */
         void timestep();

         /**
          * @brief Pure virtual method for initialising equation
          *
          * The implementation will be left open until a specialised equation is defined
          */
         virtual void init() = 0;

         /**
          * @brief Add a boundary condition to the Toroidal component equation
          *
          * \param pBC Boundary condition to add to toroidal field list
          */
         void addTorBC(SmartBC pBC);

         /**
          * @brief Add a boundary condition to the Poloidal component equation
          *
          * \param pBC Boundary condition to add to poloidal field list
          */
         void addPolBC(SmartBC pBC);

      protected:

         /**
          * @brief Timestepper object responsible for the Toroidal component
          */
         typename TimeTraits::Timestepper  mTorTStepper;

         /**
          * @brief Timestepper object responsible for the Poloidal component
          */
         typename TInfluenceTraits<TimeTraits>::Timestepper  mPolTStepper;

         /**
          * @brief Vector of BoundaryCondition for the Toroidal component
          */
         std::vector<SmartBC>    mTorBCs;

         /**
          * @brief Vector of BoundaryCondition for the Poloidal component
          */
         std::vector<SmartBC>    mPolBCs;

         /**
          * @brief Initialise the timestepper objects
          */
         void initTSteppers();

         /**
          * @brief  Compute the spectral expansion of the non linear terms
          *
          * \param rCurl Curl component of the transformed field (toroidal)
          * \param rCurlCurl Double curl component of the transformed field (poloidal)
          */
         void transformNTerms(ScalarType &rCurl, ScalarType &rCurlCurl);

         /**
          * @brief  Update value stored in OldNTerms
          */
         void updateOldNTerms();

      private:
   };

   template <typename TSimType, typename TFieldType, template <typename> class TInfluenceTraits> TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::TorPolDiffusionEquation(TFieldType &rF, typename TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::TransformType &transform, TimestepParameters &tsteps, int nBCT, int nBCP, EPMFloat a, EPMFloat b)
      : TimeEquation<TSimType, TFieldType>(rF, transform, tsteps), mTorTStepper(a, b, transform.radBasis(), tsteps, rF.trunc()), mPolTStepper(a, b, transform.radBasis(), tsteps, rF.trunc())
   {
      // Set counter to total number of boundary conditions (doesn't make the difference between toroidal or poloidal)
      this->mBCCounter = nBCT + nBCP;
   }

   template <typename TSimType, typename TFieldType, template <typename> class TInfluenceTraits> void TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::addTorBC(SmartBC pBC)
   {
      // Add toroidal BC to list
      this->mTorBCs.push_back(pBC);

      // Decrement BC counter
      --this->mBCCounter;
   }

   template <typename TSimType, typename TFieldType, template <typename> class TInfluenceTraits> void TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::addPolBC(SmartBC pBC)
   {
      // Add poloidal BC to list
      this->mPolBCs.push_back(pBC);

      // Decrement BC counter
      --this->mBCCounter;
   }

   template <typename TSimType, typename TFieldType, template <typename> class TInfluenceTraits> inline void TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::initTSteppers()
   {
      // Set Boundary condition of timestepper for Toroidal component
      for(unsigned int i=0; i < this->mTorBCs.size(); ++i)
      {
         this->mTorTStepper.addBC(this->mTorBCs.at(i));
      }

      // Initialise toroidal timestepper
      this->mTorTStepper.init();

      // Clear list of toroidal boundary conditions
      this->mTorBCs.clear();


      // Set Boundary condition of timestepper for Poloidal component
      for(unsigned int i=0; i < this->mPolBCs.size(); ++i)
      {
         this->mPolTStepper.addBC(this->mPolBCs.at(i));
      }

      // Initialise poloidal timestepper
      this->mPolTStepper.init();

      // Clear list of poloidal boundary conditions
      this->mPolBCs.clear();
   }

   template <typename TSimType, typename TFieldType, template <typename> class TInfluenceTraits> inline void TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::transformNTerms(typename TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::ScalarType &rCurl, typename TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::ScalarType &rCurlCurl)
   {
      this->mrTransform.transformRTP2TorPolNTerms(rCurl, rCurlCurl, this->mNTerms.oc().rtp());
   }

   template <typename TSimType, typename TFieldType, template <typename> class TInfluenceTraits> inline void TorPolDiffusionEquation<TSimType, TFieldType, TInfluenceTraits>::timestep()
   {
      // Timestep toroidal part
      this->mTorTStepper.timestep(this->mrX.rOc().rPerturbation().rTor(), this->mNTerms.rOc().rPerturbation().rTor());

      // Timestep poloidal part
      this->mPolTStepper.timestep(this->mrX.rOc().rPerturbation().rPol(), this->mNTerms.rOc().rPerturbation().rPol());
   }

}

#endif // TORPOLDIFFUSIONEQUATION_HPP
