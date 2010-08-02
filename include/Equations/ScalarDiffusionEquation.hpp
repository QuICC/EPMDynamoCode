/** \file ScalarDiffusionEquation.hpp
 *  \brief General scalar field diffusion equation
 */

#ifndef SCALARDIFFUSIONEQUATION_HPP
#define SCALARDIFFUSIONEQUATION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/TimeEquation.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief This class specialises the TimeEquation class to the case where the unknown
    *        is a scalar field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    */
   template <typename TSimType, typename TFieldType> class ScalarDiffusionEquation: public TimeEquation<TSimType, TFieldType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local transform type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from the simulation trait to local timestepping scheme
         typedef typename SimulationTraits<TSimType>::TimestepTraits  TimeTraits;

         /**
          * \brief Constructs the Time equation from the given general Scalar, the transform reference and the timestepping parameters
          *
          * \param rC Reference to the unknown field
          * \param transform Reference to the transform object
          * \param tsteps Timestep parameters
          * \param nBC Number of boundary conditions
          * \param a \f$a\f$ parameter for the timestepper
          * \param b \f$b\f$ parameter for the timestepper
          */
         ScalarDiffusionEquation(TFieldType &rC, TransformType &transform, TimestepParameters& tsteps, int nBC, EPMFloat a, EPMFloat b);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ScalarDiffusionEquation() {};

         /**
          * @brief Timestep the equation
          */
         void timestep();

         /**
          * @brief Pure virtual method for initialising equation
          *
          * The implementation will be left open until a specialised equation is defined
          */
         virtual void init() = 0;

         /**
          * @brief Add a boundary condition
          *
          * \param pBC Smart pointer of the boundary condition to add
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Get the maximum number of forward SSH transform packs required
          *
          * This is independent of the exact equation setup
          */
         int nFSSHPacks() const;

         /**
          * @brief Get the maximum number of forward SH transform packs required
          *
          * This is independent of the exact equation setup
          */
         int nFSHPacks() const;

         /**
          * @brief Get the maximum number of backward SSH transform packs required
          *
          * The actual number of required packs depends on exact definition of equation
          */
         virtual int nSSHBPacks() const = 0;

         /**
          * @brief Get the maximum number of backward SH transform packs required
          *
          * The actual number of required packs depends on exact definition of equation
          */
         virtual int nSHBPacks() const = 0;
         
      protected:

         /**
          * @brief Timestepper object responsible for scalar
          */
         typename TimeTraits::Timestepper  mTStepper;

         /**
          * @brief Vector of BoundaryCondition
          */
         std::vector<SmartBC>    mBCs;

         /**
          * @brief Initialise the timestepper objects
          */
         void initTStepper();

         /**
          * @brief  Compute the spectral expansion of the non linear terms
          */
         void transformNTerms();

      private:
   };

   template <typename TSimType, typename TFieldType> inline int ScalarDiffusionEquation<TSimType,TFieldType>::nFSSHPacks() const
   {
      return 1;
   }

   template <typename TSimType, typename TFieldType> inline int ScalarDiffusionEquation<TSimType,TFieldType>::nFSHPacks() const
   {
      return 1;
   }

   template <typename TSimType, typename TFieldType> ScalarDiffusionEquation<TSimType,TFieldType>::ScalarDiffusionEquation(TFieldType &rC, typename ScalarDiffusionEquation<TSimType, TFieldType>::TransformType &transform, TimestepParameters& tsteps, int nBC, EPMFloat a, EPMFloat b)
      : TimeEquation<TSimType, TFieldType>(rC, transform, tsteps), mTStepper(a, b, transform.radBasis(), tsteps, rC.oc().trunc(), true)
   {
      // Set the boundary condition counter
      this->mBCCounter = nBC;
   }

   template <typename TSimType, typename TFieldType> void ScalarDiffusionEquation<TSimType, TFieldType>::addBC(SmartBC pBC)
   {
      // Add boundary condition to list
      this->mBCs.push_back(pBC);

      // Decrement the BC counter
      --this->mBCCounter;
   }

   template <typename TSimType, typename TFieldType> inline void ScalarDiffusionEquation<TSimType, TFieldType>::initTStepper()
   {
      // Set Boundary condition of timestepper
      for(unsigned int i=0; i < this->mBCs.size(); ++i)
      {
         this->mTStepper.addBC(this->mBCs.at(i));
      }
      
      // Initialise the timestepper
      this->mTStepper.init();

      // Clear the list of boundary conditions
      this->mBCs.clear();
   }

   template <typename TSimType, typename TFieldType> inline void ScalarDiffusionEquation<TSimType, TFieldType>::transformNTerms()
   {
      // Transform the real space non linear terms to spectral space
      this->mrTransform.transformRTP2Spec(this->mNTerms.rOc().rPerturbation(), this->mNTerms.oc().rtp());
   }

   template <typename TSimType, typename TFieldType> inline void ScalarDiffusionEquation<TSimType, TFieldType>::timestep()
   {
      // Timestep equation
      this->mTStepper.timestep(this->mrX.rOc().rPerturbation(), this->mNTerms.rOc().rPerturbation());
   }

}

#endif // SCALARDIFFUSIONEQUATION_HPP
