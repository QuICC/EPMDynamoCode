/** \file NavierStokesLinRotating.hpp
 *  \brief Implementation of the Navier-Stokes equation for rotating convection
 */

#ifndef NAVIERSTOKESLINROTATING_HPP
#define NAVIERSTOKESLINROTATING_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/TorPolDiffusionEquation.hpp"
#include "Equations/NavierStokes/NavierStokesThermal.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with rotating convection
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesLinRotating : public NavierStokesThermal<TSimType, TSimTraits>  {
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
         NavierStokesLinRotating(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesLinRotating() {};

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

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesLinRotating<TSimType, TSimTraits>::NavierStokesLinRotating(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, typename NavierStokesLinRotating<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesLinRotating<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesThermal<TSimType, TSimTraits>(rV, rC, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesLinRotating<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesLinRotating<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<0>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESLINROTATING_HPP
