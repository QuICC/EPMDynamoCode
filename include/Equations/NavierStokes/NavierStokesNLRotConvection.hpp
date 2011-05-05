/** \file NavierStokesNLRotConvection.hpp
 *  \brief Implementation of a rotating convection equation
 */

#ifndef NAVIERSTOKESNLROTCONVECTION_HPP
#define NAVIERSTOKESNLROTCONVECTION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Config/SimulationConfig.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesConvection.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with rotating convection
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesNLRotConvection : public NavierStokesConvection<TSimType, TSimTraits>
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
          * @param rC Codensity scalar (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesNLRotConvection(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesNLRotConvection() {};

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

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesNLRotConvection<TSimType, TSimTraits>::NavierStokesNLRotConvection(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, typename NavierStokesNLRotConvection<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesNLRotConvection<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesConvection<TSimType, TSimTraits>(rV, rC, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesNLRotConvection<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesNLRotConvection<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());

      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESNLROTCONVECTION_HPP
