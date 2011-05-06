/** \file NavierStokesMHD.hpp
 *  \brief Implementation of the MHD Navier-Stokes equation
 */

#ifndef NAVIERSTOKESMHD_HPP
#define NAVIERSTOKESMHD_HPP

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
    * @brief This class implements the full MHD Navierr-Stokes equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class NavierStokesMHD : public NavierStokesConvection<TSimTraits>
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
          * @param rB Magnetic field (stored as reference)
          * @param rC Codensity scalar (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesMHD(typename TSimTraits::VelType &rV, typename TSimTraits::MagType &rB, typename TSimTraits::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesMHD() {};

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
          * @brief Const Reference variable to the magnetic field
          */
         typename TSimTraits::MagType&  mrB;

      private:
   };

   template <typename TSimTraits> NavierStokesMHD<TSimTraits>::NavierStokesMHD(typename TSimTraits::VelType &rV, typename TSimTraits::MagType &rB, typename TSimTraits::CodType &rC, typename NavierStokesMHD<TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesMHD<TSimTraits>::EquationParametersType &params)
      : NavierStokesConvection<TSimTraits>(rV, rC, transform, tsteps, params), mrB(rB)
   {
   }

   template <typename TSimTraits> void NavierStokesMHD<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);

      // Update real space values of magnetic field
      this->mrB.rOc().transform(step);

      // Update real space values of the curl of the magnetic field
      this->mrB.rOc().curlTransform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimTraits> void NavierStokesMHD<TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());

      // Compute \f$(\nabla \times B)\times B\f$
      this->mrB.oc().curl().template cross<1>(this->mNTerms.rOc().rRTP(), this->mrB.oc().rtp(), this->mrParams.nsLorentz());

      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESMHD_HPP
