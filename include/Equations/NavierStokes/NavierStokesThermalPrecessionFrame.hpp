/** \file NavierStokesThermalPrecessionFrame.hpp
 *  General representation of the Navier-Stokes diffusion equation with non-linear and precession term in the precession frame
 */

#ifndef NAVIERSTOKESTHERMALPRECESSIONFRAME_HPP
#define NAVIERSTOKESTHERMALPRECESSIONFRAME_HPP

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
#include "General/MathConstants.hpp"
#include "Equations/NavierStokes/NavierStokesDiffusion.hpp"
#include "IO/ASCII/ConfigurationFile.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes diffusion equation with non-linear and precession term in the precession frame
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class NavierStokesThermalPrecessionFrame : public NavierStokesDiffusion<TSimTraits>
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
          * @param rC Codensity scakar  (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesThermalPrecessionFrame(typename TSimTraits::VelType &rV, typename TSimTraits::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesThermalPrecessionFrame() {};

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
         * @brief Const reference vatiable to the codensity scalar
         *
         */         
       
         typename TSimTraits::CodType& mrC;

         /**
          * @brief Const Reference variable to the magnetic field
          */
         //typename TSimTraits::MagType&  mrB;


      private:
         /**
          * @brief Setup precession forcing
          */
         void setupPrecession();

         /**
          * @brief Precession rotation rate
          */
         EPMFloat mOmega;

         /**
          * @brief Precession tilt angle cosinus
          */
         EPMFloat mCosAlpha;

         /**
          * @brief Precession tilt angle sinus
          */
         EPMFloat mSinAlpha;
   };

   template <typename TSimTraits> NavierStokesThermalPrecessionFrame<TSimTraits>::NavierStokesThermalPrecessionFrame(typename TSimTraits::VelType &rV,  typename TSimTraits::CodType &rC, typename NavierStokesThermalPrecessionFrame<TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesThermalPrecessionFrame<TSimTraits>::EquationParametersType &params)
      : NavierStokesDiffusion<TSimTraits>(rV, transform, tsteps, params), mrC(rC), mOmega(0.0), mCosAlpha(0.0), mSinAlpha(0.0)
   {
      // Setup precession forcing
      this->setupPrecession();
   }

   template <typename TSimTraits> void NavierStokesThermalPrecessionFrame<TSimTraits>::setupPrecession()
   {
      std::vector<std::string>   integers;
      std::vector<std::string>   floats;

      floats.push_back("omega");
      floats.push_back("alpha");

      ConfigurationFile cfg("precession", integers, floats);

      // Initialise configuration file
      cfg.init();

      // Read data
      cfg.read();

      // Finalise configuration file reader
      cfg.finalise();

      // Print information
      cfg.printInfo();

      // Store values for configuration file
      this->mOmega = cfg.floats()(0);
      // read angle and convert to radians
      this->mCosAlpha = std::cos(cfg.floats()(1)*(MathConstants::PI/180.));
      this->mSinAlpha = std::sin(cfg.floats()(1)*(MathConstants::PI/180.));
   }

   template <typename TSimTraits> void NavierStokesThermalPrecessionFrame<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);

      // Update real space values of magnetic field
      //this->mrB.rOc().transform(step);

      // Update real space values of the curl of the magnetic field
      //this->mrB.rOc().curlTransform(step);
   }

   template <typename TSimTraits> void NavierStokesThermalPrecessionFrame<TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.nsAdvection());

      // Compute \f$C \vec{r}\f$
      // Carefull with radVec<0> sets the value instead of adding
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());

      // Compute \f$(\nabla \times B)\times B\f$
      //this->mrB.oc().curl().template cross<1>(this->mNTerms.rOc().rRTP(), this->mrB.oc().rtp(), this->mrParams.nsLorentz());

      // Compute \f$\Omega_p\times\vec{u}\f$
      this->mrX.oc().rtp().template precessionFrame<-1>(this->mNTerms.rOc().rRTP(), this->mOmega, this->mCosAlpha, this->mSinAlpha, this->mrParams.nsCoriolis());

      // Compute poincare force \f$\Omega\times\hat{z}\vec{r}\f$
      //RTPOperators::subPoincare(this->mNTerms.rOc().rRTP(), this->mOmega, this->mCosAlpha, this->mSinAlpha, this->mrTStepParams.time());
      //Leo: Thermal convection with diffusion time scale
      //RTPOperators::subPoincare(this->mNTerms.rOc().rRTP(), this->mOmega, this->mCosAlpha, this->mSinAlpha, this->mrTStepParams.time(), 1./this->mrParams.E());
   }

}

#endif // NAVIERSTOKESTHERMALPRECESSIONFRAME_HPP
