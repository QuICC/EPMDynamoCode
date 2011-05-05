/** \file NavierStokesDiffusion.hpp
 *  \brief Implementation of a velocity diffusion equation
 */

#ifndef NAVIERSTOKESDIFFUSION_HPP
#define NAVIERSTOKESDIFFUSION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Config/SimulationConfig.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesBase.hpp"
#include "Timestepping/Traits/InfluenceTTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with only diffusion
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesDiffusion : public NavierStokesBase<TSimType, TSimTraits, InfluenceTTraits>
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
         NavierStokesDiffusion(typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

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

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesDiffusion<TSimType, TSimTraits>::NavierStokesDiffusion(typename TSimTraits<TSimType>::VelType &rV, typename NavierStokesDiffusion<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesDiffusion<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesBase<TSimType, TSimTraits, InfluenceTTraits>(rV, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesDiffusion<TSimType, TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesDiffusion<TSimType, TSimTraits>::updateRHS()
   {
      // Set the RTP non linear terms to zero
      this->mNTerms.rOc().rRTP().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesDiffusion<TSimType, TSimTraits>::transformRHS(const int step)
   {
      // Transform non linear terms to spectral space from mNTerms values
      this->transformNTerms(this->mNTerms.rOc().rPerturbation().rTor(), this->mNTerms.rOc().rPerturbation().rPol());
   }

}

#endif // NAVIERSTOKESDIFFUSION_HPP
