/** \file NavierStokesThermal.hpp
 *  \brief Implementation of the Navier-Stokes equation for thermal convection
 */

#ifndef NAVIERSTOKESTHERMAL_HPP
#define NAVIERSTOKESTHERMAL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/TorPolTimeEquation.hpp"
#include "Equations/NavierStokes/NavierStokesRotating.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation with rotation and thermal convection
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesThermal : public NavierStokesRotating<TSimType, TSimTraits>  {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSim>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param rV Magnetic field (stored as reference)
          * @param rC Codensity scalar (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesThermal(typename TSimTraits<TSim>::VelType &rV, typename TSimTraits<TSim>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParameters &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesThermal();

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
         typename TSimTraits<TSim>::CodType&  mrC;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesThermal<TSimType, TSimTraits>::NavierStokesThermal(typename TSimTraits<TSim>::VelType &rV, typename TSimTraits<TSim>::CodType &rC, typename NavierStokesThermal<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, EquationParameters &params)
      : NavierStokesRotating<TSimType, TSimTraits>(rV, transform, tsteps, params), mrC(rC)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesThermal<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesThermal<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.Ro());

      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.q()*this->mrParams.Ra());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP());
   }

}

#endif // NAVIERSTOKESTHERMAL_HPP
