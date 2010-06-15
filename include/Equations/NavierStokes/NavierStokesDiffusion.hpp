/** \file NavierStokesDiffusion.hpp
 *  \brief Implementation of the velocity diffusion equation
 */

#ifndef NAVIERSTOKESDIFFUSION_HPP
#define NAVIERSTOKESDIFFUSION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesBase.hpp"
#include "Equations/EquationParameters.hpp"
#include "Timestepping/WithoutPoissonTStep.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with only diffusion
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    * \tparam TTStepPol NEED TO CHANGE THIS
    *
    * \bug Need to change the implementation of the poisson solver
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesDiffusion : public NavierStokesBase<TSimType, TSimTraits, WithoutPoissonTStep> 
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSim>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesDiffusion(typename TSimTraits<TSim>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParameters &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesDiffusion() {};

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
         
      protected:

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesDiffusion<TSimType, TSimTraits>::NavierStokesDiffusion(typename TSimTraits<TSim>::VelType &rV, typename NavierStokesDiffusion<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, EquationParameters &params)
      : NavierStokesBase<TSimType, TSimTraits, WithoutPoissonTStep>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesDiffusion<TSimType, TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesDiffusion<TSimType, TSimTraits>::updateRHS()
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesDiffusion<TSimType, TSimTraits>::transformRHS(const int step)
   {
      if(step == 0)
      {
         // Set the non linear terms to zero
         int nL = this->mNTerms.oc().perturbation().nL();
         const int l0 = this->mNTerms.oc().perturbation().tor().minL();

         for(int l = l0; l < nL; ++l)
         {
            this->mNTerms.rOc().rPerturbation().rTor().rLShell(l).setConstant(0.0);
            this->mNTerms.rOc().rPerturbation().rPol().rLShell(l).setConstant(0.0);
         }
      }
   }

}

#endif // NAVIERSTOKESDIFFUSION_HPP
