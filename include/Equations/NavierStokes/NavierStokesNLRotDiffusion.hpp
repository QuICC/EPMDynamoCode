/** \file NavierStokesNLRotDiffusion.hpp
 *  General representation of the Navier-Stokes diffusion equation with non-linear and coriolis term
 */

#ifndef NAVIERSTOKESNLROTDIFFUSION_HPP
#define NAVIERSTOKESNLROTDIFFUSION_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesDiffusion.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes diffusion equation with non-linear and coriolis term
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class NavierStokesNLRotDiffusion : public NavierStokesDiffusion<TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for the EquationParameters type
         typedef SimulationConfig::EquationParametersType EquationParametersType;

         /**
          * @brief Constructor
          *
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesNLRotDiffusion(typename TSimTraits::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

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

   template <typename TSimTraits> NavierStokesNLRotDiffusion<TSimTraits>::NavierStokesNLRotDiffusion(typename TSimTraits::VelType &rV, typename NavierStokesNLRotDiffusion<TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesNLRotDiffusion<TSimTraits>::EquationParametersType &params)
      : NavierStokesDiffusion<TSimTraits>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimTraits> void NavierStokesNLRotDiffusion<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);
   }

   template <typename TSimTraits> void NavierStokesNLRotDiffusion<TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESNLROTDIFFUSION_HPP
