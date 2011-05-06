/** \file NavierStokesNLConvection.hpp
 *  \brief Implementation of a rotating convection equation
 */

#ifndef NAVIERSTOKESNLCONVECTION_HPP
#define NAVIERSTOKESNLCONVECTION_HPP

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
#include "Equations/NavierStokes/NavierStokesConvection.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with rotating convection
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class NavierStokesNLConvection : public NavierStokesConvection<TSimTraits>
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
         NavierStokesNLConvection(typename TSimTraits::VelType &rV, typename TSimTraits::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesNLConvection() {};

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

   template <typename TSimTraits> NavierStokesNLConvection<TSimTraits>::NavierStokesNLConvection(typename TSimTraits::VelType &rV, typename TSimTraits::CodType &rC, typename NavierStokesNLConvection<TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesNLConvection<TSimTraits>::EquationParametersType &params)
      : NavierStokesConvection<TSimTraits>(rV, rC, transform, tsteps, params)
   {
   }

   template <typename TSimTraits> void NavierStokesNLConvection<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimTraits> void NavierStokesNLConvection<TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());

      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());
   }

}

#endif // NAVIERSTOKESNLCONVECTION_HPP
