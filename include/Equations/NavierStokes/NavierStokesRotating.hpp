/** \file NavierStokesRotating.hpp
 *  \brief Implementation of a rotating diffusion equation
 */

#ifndef NAVIERSTOKESROTATING_HPP
#define NAVIERSTOKESROTATING_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesBase.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "Timestepping/InfluenceTTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with rotation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesRotating : public NavierStokesBase<TSimType, TSimTraits, InfluenceTTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param rV Magnetic field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesRotating(typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParameters &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesRotating() {};

         /**
          * @brief Update RTP values of the equation
          *
          * \param step Current step in a multistep transform
          */
         void updateRTP(const int step);

         /**
          * @brief Update RHS of the equation
          */
         void updateRHS();

         /**
          * @brief Transform RHS of the equation
          *
          * \param step Current step in a multistep transform
          */
         void transformRHS(const int step);

         /**
          * @brief Initialise the equation
          */
         void init();
         
      protected:

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesRotating<TSimType, TSimTraits>::NavierStokesRotating(typename TSimTraits<TSimType>::VelType &rV, typename NavierStokesRotating<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, EquationParameters &params)
      : NavierStokesBase<TSimType, TSimTraits, InfluenceTTraits>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesRotating<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesRotating<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.Ro());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP());
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesRotating<TSimType, TSimTraits>::transformRHS(const int step)
   {
      // Transform non linear terms to spectral space from mNTerms values
      this->transformNTerms(this->mNTerms.rOc().rPerturbation().rTor(), this->mNTerms.rOc().rPerturbation().rPol());
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesRotating<TSimType, TSimTraits>::init()
   {
      // Set influence matrix equation boundary condition
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mrTransform.radBasis()));
      this->mPolTStepper.addInfluenceBC(pZeroBC);

      // Check that the right number of BCs have been provided
      if(this->hasAllBCs())
      {
         // Init the timesteppers (they depend on the boundary conditions)
         this->initTSteppers();
      }
      else
      {
         throw EPMException("NavierStokesRotating::init", "Tried to initialise with wrong number of BCs");
      }
   }

}

#endif // NAVIERSTOKESROTATING_HPP
