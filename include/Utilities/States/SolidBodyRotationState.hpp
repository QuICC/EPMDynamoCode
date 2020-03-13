/** \file SolidBodyRotationState.hpp
 *  \brief Initial state definition function for a random perturbation state
 */

#ifndef SOLIDBODYROTATIONSTATE_HPP
#define SOLIDBODYROTATIONSTATE_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits of random perturbation state generator
    */
   class SolidBodyRotationTraits
   {
      public:
         /// Requires RTP Codensity computations
         static const bool UseRTPCodensity = false;

         /// Requires spectral Codensity computations
         static const bool UseSpecCodensity = false;

         /// Requires Codensity gradient computations
         static const bool UseCodensityGrad = false;

         /// Requires RTP Magnetic computations
         static const bool UseRTPMagnetic = false;

         /// Requires spectral Magnetic computations
         static const bool UseSpecMagnetic = false;

         /// Requires Magnetic curl computations
         static const bool UseMagneticCurl = false;

         /// Requires RTP Velocity computations
         static const bool UseRTPVelocity = true;

         /// Requires Velocity computations
         static const bool UseSpecVelocity = false;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;
   };

   /**
    * \brief Initial state definition function for a random perturbation state
    *
    * \epmBug State generation is not aware of parallelisation
    */
   template <typename TGenTraits> class SolidBodyRotationState
   {
      public:
         /// Typdef for the StateTraits type
         typedef SolidBodyRotationTraits  StateTraits;

         /**
          * @brief Perturbation amplitude
          */
         static const EPMFloat PERTURBATION_AMPLITUDE;

         /**
          * @brief Ratio of L spectrum to be perturbed
          */
         static const int PERTURBATION_LRATIO;

         /**
          * @brief Ratio of radial spectrum to be perturbed
          */
         static const int PERTURBATION_NRATIO;

         /// Typedef for the codensity type
         typedef typename TGenTraits::CodType  Codensity;

         /// Typedef for the magnetic type
         typedef typename TGenTraits::MagType  Magnetic;

         /// Typedef for the velocity type
         typedef typename TGenTraits::VelType  Velocity;

         /**
          * @brief Set RTP value for the codensity scalar
          */
         static void setRTPCodensity(Codensity &codC);

         /**
          * @brief Set RTP value for the magnetic field
          */
         static void setRTPMagnetic(Magnetic &magB);

         /**
          * @brief Set RTP value for the velocity field
          */
         static void setRTPVelocity(Velocity &velV);

         /**
          * @brief Set spectral coefficients for the codensity scalar
          */
         static void setSpecCodensity(Codensity &codC);

         /**
          * @brief Set spectral coefficients for the velocity field
          */
         static void setSpecMagnetic(Magnetic &magB);

         /**
          * @brief Set spectral coefficients for the magnetic field
          */
         static void setSpecVelocity(Velocity &velV);

      private:
         /**
          * @brief Private constructor
          */
         SolidBodyRotationState();

         /**
          * @brief Destructor
          */
         virtual ~SolidBodyRotationState() {};
   };

   template <typename TGenTraits> const EPMFloat SolidBodyRotationState<TGenTraits>::PERTURBATION_AMPLITUDE = 1.0e-10;

   template <typename TGenTraits> const int SolidBodyRotationState<TGenTraits>::PERTURBATION_LRATIO = 2;

   template <typename TGenTraits> const int SolidBodyRotationState<TGenTraits>::PERTURBATION_NRATIO = 4;

   template <typename TGenTraits> void SolidBodyRotationState<TGenTraits>::setRTPCodensity(typename SolidBodyRotationState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void SolidBodyRotationState<TGenTraits>::setRTPMagnetic(typename SolidBodyRotationState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      Array sPh = pTrunc->sim()->hoz()->phGrid().array().sin();
      Array cPh = pTrunc->sim()->hoz()->phGrid().array().cos();

      EPMFloat xAmp = 1.0;
      EPMFloat yAmp = 1.0;
      EPMFloat zAmp = 1.0;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            // Create solid body rotation around x axis
            magB.rOc().rRTP().rTheta().rShell(n).col(th) += -xAmp*pTrunc->local()->rtp()->radGrid(n) * sPh;
            magB.rOc().rRTP().rPhi().rShell(n).col(th) += -xAmp*pTrunc->local()->rtp()->radGrid(n)*pTrunc->local()->rtp()->cTh(th, n) * cPh;

            // Add solid body rotation around y axis
            magB.rOc().rRTP().rTheta().rShell(n).col(th) += yAmp*pTrunc->local()->rtp()->radGrid(n) * cPh;
            magB.rOc().rRTP().rPhi().rShell(n).col(th) += -yAmp*pTrunc->local()->rtp()->radGrid(n)*pTrunc->local()->rtp()->cTh(th, n) * sPh;

            // Add solid body rotation around z axis
            magB.rOc().rRTP().rPhi().rShell(n).col(th).array() += zAmp*pTrunc->local()->rtp()->radGrid(n)*pTrunc->local()->rtp()->sTh(th, n);
         }
      }
   }

   template <typename TGenTraits> void SolidBodyRotationState<TGenTraits>::setRTPVelocity(typename SolidBodyRotationState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      Array sPh = pTrunc->sim()->hoz()->phGrid().array().sin();
      Array cPh = pTrunc->sim()->hoz()->phGrid().array().cos();

      EPMFloat xAmp = 1.0;
      EPMFloat yAmp = 1.0;
      EPMFloat zAmp = 1.0;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            // Create solid body rotation around x axis
            velV.rOc().rRTP().rTheta().rShell(n).col(th) += -xAmp*pTrunc->local()->rtp()->radGrid(n) * sPh;
            velV.rOc().rRTP().rPhi().rShell(n).col(th) += -xAmp*pTrunc->local()->rtp()->radGrid(n)*pTrunc->local()->rtp()->cTh(th, n) * cPh;

            // Add solid body rotation around y axis
            velV.rOc().rRTP().rTheta().rShell(n).col(th) += yAmp*pTrunc->local()->rtp()->radGrid(n) * cPh;
            velV.rOc().rRTP().rPhi().rShell(n).col(th) += -yAmp*pTrunc->local()->rtp()->radGrid(n)*pTrunc->local()->rtp()->cTh(th, n) * sPh;

            // Add solid body rotation around z axis
            velV.rOc().rRTP().rPhi().rShell(n).col(th).array() += zAmp*pTrunc->local()->rtp()->radGrid(n)*pTrunc->local()->rtp()->sTh(th, n);
         }
      }
   }

   template <typename TGenTraits> void SolidBodyRotationState<TGenTraits>::setSpecCodensity(typename SolidBodyRotationState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      // Set some perturbation random energy
      for(int l=0; l < pTrunc->local()->spec()->nL()/SolidBodyRotationState<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         codC.rOc().rPerturbation().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/SolidBodyRotationState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         codC.rOc().rPerturbation().rLShell(l) *= EPMComplex(SolidBodyRotationState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            codC.rOc().rPerturbation().rLShell(l).col(0)(n).imag(0.0);
         }
      }

      // Create basic state
      //codC.rOc().rPerturbation().rLShell(0)(0,0) += 1.0/4.0;
      codC.rOc().rPerturbation().rLShell(0)(0,0) += 0.313329;
      //codC.rOc().rPerturbation().rLShell(0)(1,0) += -1.0/2.0;
      codC.rOc().rPerturbation().rLShell(0)(1,0) += -0.221557;
   }

   template <typename TGenTraits> void SolidBodyRotationState<TGenTraits>::setSpecMagnetic(typename SolidBodyRotationState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/SolidBodyRotationState<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         // Set some perturbation random energy in Toroidal component
         magB.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/SolidBodyRotationState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(SolidBodyRotationState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Set some perturbation random energy in Poloidal component
         magB.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/SolidBodyRotationState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(SolidBodyRotationState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag(0.0);
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag(0.0);
         }
      }
   }

   template <typename TGenTraits> void SolidBodyRotationState<TGenTraits>::setSpecVelocity(typename SolidBodyRotationState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/SolidBodyRotationState<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         // Set some perturbation random energy in Toroidal component
         velV.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/SolidBodyRotationState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         velV.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(SolidBodyRotationState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Set some perturbation random energy in Poloidal component
         velV.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/SolidBodyRotationState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         velV.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(SolidBodyRotationState<TGenTraits>::PERTURBATION_AMPLITUDE, 0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            velV.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag(0.0);
            velV.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag(0.0);
         }
      }

   }

   template <typename TGenTraits> SolidBodyRotationState<TGenTraits>::SolidBodyRotationState()
   {
   }
}

#endif // SOLIDBODYROTATIONSTATE_HPP
