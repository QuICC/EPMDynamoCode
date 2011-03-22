/** \file NavierStokesRotDiffusion.hpp
 *  General representation of the Navier-Stokes diffusion equation with coriolis term
 */

#ifndef NAVIERSTOKESROTDIFFUSION_HPP
#define NAVIERSTOKESROTDIFFUSION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesDiffusion.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes diffusion equation with coriolis term
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesRotDiffusion : public NavierStokesDiffusion<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for the EquationParameters type
         typedef typename SimulationTraits<TSimType>::EquationParametersType EquationParametersType;

         /**
          * @brief Constructor
          *
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesRotDiffusion(typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesRotDiffusion() {};

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
         
      protected:

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesRotDiffusion<TSimType, TSimTraits>::NavierStokesRotDiffusion(typename TSimTraits<TSimType>::VelType &rV, typename NavierStokesRotDiffusion<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesRotDiffusion<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesDiffusion<TSimType, TSimTraits>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesRotDiffusion<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesRotDiffusion<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<0>(this->mNTerms.rOc().rRTP(), -this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESROTDIFFUSION_HPP
