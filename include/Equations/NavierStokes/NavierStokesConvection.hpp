/** \file NavierStokesConvection.hpp
 *  \brief Implementation of a rotating convection equation
 */

#ifndef NAVIERSTOKESCONVECTION_HPP
#define NAVIERSTOKESCONVECTION_HPP

// Configuration includes
//

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
    * @brief General representation of the Navier-Stokes equation with rotating convection
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class NavierStokesConvection : public NavierStokesDiffusion<TSimTraits>
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
          * @param rC Codensity scalar (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesConvection(typename TSimTraits::VelType &rV, typename TSimTraits::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesConvection() {};

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

         /**
          * @brief Const Reference variable to the codensity scalar
          */
         typename TSimTraits::CodType&  mrC;

      private:
   };

   template <typename TSimTraits> NavierStokesConvection<TSimTraits>::NavierStokesConvection(typename TSimTraits::VelType &rV, typename TSimTraits::CodType &rC, typename NavierStokesConvection<TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesConvection<TSimTraits>::EquationParametersType &params)
      : NavierStokesDiffusion<TSimTraits>(rV, transform, tsteps, params), mrC(rC)
   {
   }

   template <typename TSimTraits> void NavierStokesConvection<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimTraits> void NavierStokesConvection<TSimTraits>::updateRHS()
   {
      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<0>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());
   }

}

#endif // NAVIERSTOKESCONVECTION_HPP
