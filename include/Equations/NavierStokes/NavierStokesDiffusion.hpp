/** \file NavierStokesDiffusion.hpp
 *  \brief Implementation of a velocity diffusion equation
 */

#ifndef NAVIERSTOKESDIFFUSION_HPP
#define NAVIERSTOKESDIFFUSION_HPP

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
#include "Equations/NavierStokes/NavierStokesBase.hpp"
#include "Timestepping/Traits/InfluenceTTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with only diffusion
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class NavierStokesDiffusion : public NavierStokesBase<TSimTraits, InfluenceTTraits>
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
         NavierStokesDiffusion(typename TSimTraits::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

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

   template <typename TSimTraits> NavierStokesDiffusion<TSimTraits>::NavierStokesDiffusion(typename TSimTraits::VelType &rV, typename NavierStokesDiffusion<TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesDiffusion<TSimTraits>::EquationParametersType &params)
      : NavierStokesBase<TSimTraits, InfluenceTTraits>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimTraits> void NavierStokesDiffusion<TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimTraits> void NavierStokesDiffusion<TSimTraits>::updateRHS()
   {
      // Set the RTP non linear terms to zero
      this->mNTerms.rOc().rRTP().initialiseZeros();
   }

   template <typename TSimTraits> void NavierStokesDiffusion<TSimTraits>::transformRHS(const int step)
   {
      // Transform non linear terms to spectral space from mNTerms values
      this->transformNTerms(this->mNTerms.rOc().rPerturbation().rTor(), this->mNTerms.rOc().rPerturbation().rPol());
   }

}

#endif // NAVIERSTOKESDIFFUSION_HPP
