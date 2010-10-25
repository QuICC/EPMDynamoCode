/** \file NavierStokesLinThermal.hpp
 *  \brief Implementation of the Navier-Stokes equation for thermal convection
 */

#ifndef NAVIERSTOKESLINTHERMAL_HPP
#define NAVIERSTOKESLINTHERMAL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/TorPolDiffusionEquation.hpp"
#include "Equations/NavierStokes/NavierStokesRotating.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with rotation and thermal convection
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesLinThermal : public NavierStokesRotating<TSimType, TSimTraits>  {
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
          * @param rV Magnetic field (stored as reference)
          * @param rC Codensity scalar (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesLinThermal(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesLinThermal() {};

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
         typename TSimTraits<TSimType>::CodType&  mrC;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesLinThermal<TSimType, TSimTraits>::NavierStokesLinThermal(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, typename NavierStokesLinThermal<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesLinThermal<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesRotating<TSimType, TSimTraits>(rV, transform, tsteps, params), mrC(rC)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesLinThermal<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesLinThermal<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<0>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESLINTHERMAL_HPP
