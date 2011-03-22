/** \file NavierStokesNLRotDiffusion.hpp
 *  General representation of the Navier-Stokes diffusion equation with non-linear and coriolis term
 */

#ifndef NAVIERSTOKESNLROTDIFFUSION_HPP
#define NAVIERSTOKESNLROTDIFFUSION_HPP

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
    * @brief General representation of the Navier-Stokes diffusion equation with non-linear and coriolis term
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesNLRotDiffusion : public NavierStokesDiffusion<TSimType, TSimTraits>
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
         NavierStokesNLRotDiffusion(typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesNLRotDiffusion() {};

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

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesNLRotDiffusion<TSimType, TSimTraits>::NavierStokesNLRotDiffusion(typename TSimTraits<TSimType>::VelType &rV, typename NavierStokesNLRotDiffusion<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesNLRotDiffusion<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesDiffusion<TSimType, TSimTraits>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesNLRotDiffusion<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesNLRotDiffusion<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESNLROTDIFFUSION_HPP
