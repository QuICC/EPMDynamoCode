/** \file NavierStokesNLDiffusion.hpp
 *  General representation of the Navier-Stokes diffusion equation with non-linear term
 */

#ifndef NAVIERSTOKESNLDIFFUSION_HPP
#define NAVIERSTOKESNLDIFFUSION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Config/SimulationConfig.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesDiffusion.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes diffusion equation with non-linear term
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesNLDiffusion : public NavierStokesDiffusion<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

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
         NavierStokesNLDiffusion(typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesNLDiffusion() {};

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

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesNLDiffusion<TSimType, TSimTraits>::NavierStokesNLDiffusion(typename TSimTraits<TSimType>::VelType &rV, typename NavierStokesNLDiffusion<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesNLDiffusion<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesDiffusion<TSimType, TSimTraits>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesNLDiffusion<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesNLDiffusion<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());
   }

}

#endif // NAVIERSTOKESNLDIFFUSION_HPP
