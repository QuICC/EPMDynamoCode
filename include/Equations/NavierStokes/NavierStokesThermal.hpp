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
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/TorPolDiffusionEquation.hpp"
#include "Equations/NavierStokes/NavierStokesBase.hpp"
#include "Timestepping/Traits/InfluenceTTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation thermal convection
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesThermal : public NavierStokesBase<TSimType, TSimTraits, InfluenceTTraits>  {
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
         NavierStokesThermal(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesThermal() {};

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

         /**
          * @brief Const Reference variable to the codensity scalar
          */
         typename TSimTraits<TSimType>::CodType&  mrC;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesThermal<TSimType, TSimTraits>::NavierStokesThermal(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::CodType &rC, typename NavierStokesThermal<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesThermal<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesBase<TSimType, TSimTraits, InfluenceTTraits>(rV, transform, tsteps, params), mrC(rC)
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
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());

      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesThermal<TSimType, TSimTraits>::transformRHS(const int step)
   {
      // Transform non linear terms to spectral space from mNTerms values
      this->transformNTerms(this->mNTerms.rOc().rPerturbation().rTor(), this->mNTerms.rOc().rPerturbation().rPol());
   }

}

#endif // NAVIERSTOKESTHERMAL_HPP
