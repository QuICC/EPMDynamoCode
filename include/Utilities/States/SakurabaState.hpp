/** \file SakurabaState.hpp
 *  \brief Initial state definition function for Sakuraba's paper
 */

#ifndef SAKURABASTATE_HPP
#define SAKURABASTATE_HPP

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
    * @brief Traits of Sakuraba state generator
    */
   class SakurabaTraits
   {
      public:
         /// Requires RTP Codensity computations
         static const bool UseRTPCodensity = false;

         /// Requires spectral Codensity computations
         static const bool UseSpecCodensity = true;

         /// Requires Codensity gradient computations
         static const bool UseCodensityGrad = false;

         /// Requires RTP Magnetic computations
         static const bool UseRTPMagnetic = false;

         /// Requires spectral Magnetic computations
         static const bool UseSpecMagnetic = true;

         /// Requires Magnetic curl computations
         static const bool UseMagneticCurl = false;

         /// Requires RTP Velocity computations
         static const bool UseRTPVelocity = false;

         /// Requires Velocity computations
         static const bool UseSpecVelocity = true;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;
   };

   /**
    * \brief Initial state definition function for Sakuraba's paper
    */
   template <typename TSimTraits> class SakurabaState
   {
      public:
         /// Typdef for the StateTraits type
         typedef SakurabaTraits  StateTraits;

         /**
          * @brief Perturb the magnetic field?
          */
         static const bool PERTURB_MAGNETIC = true;

         /**
          * @brief Perturbation amplitude
          */
         static const EPMFloat PERTURBATION_AMPLITUDE;

         /// Typedef for the codensity type
         typedef typename TSimTraits::CodType  Codensity;

         /// Typedef for the magnetic type
         typedef typename TSimTraits::MagType  Magnetic;

         /// Typedef for the velocity type
         typedef typename TSimTraits::VelType  Velocity;

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
         SakurabaState();

         /**
          * @brief Destructor
          */
         virtual ~SakurabaState() {};
   };

   template <typename TSimTraits> const EPMFloat SakurabaState<TSimTraits>::PERTURBATION_AMPLITUDE = 1e-9;

   template <typename TSimTraits> void SakurabaState<TSimTraits>::setRTPCodensity(typename SakurabaState<TSimTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TSimTraits> void SakurabaState<TSimTraits>::setRTPMagnetic(typename SakurabaState<TSimTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
      }
   }

   template <typename TSimTraits> void SakurabaState<TSimTraits>::setRTPVelocity(typename SakurabaState<TSimTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
      }
   }

   template <typename TSimTraits> void SakurabaState<TSimTraits>::setSpecCodensity(typename SakurabaState<TSimTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      // Set some perturbation random energy
      for(int l=0; l < pTrunc->local()->spec()->nL()/2; ++l)
      {
         codC.rOc().rPerturbation().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l)/2, 1)).setRandom();
         codC.rOc().rPerturbation().rLShell(l) *= EPMComplex(SakurabaState<TSimTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            codC.rOc().rPerturbation().rLShell(l).col(0)(n).imag(0.0);
         }
      }

      // Create basic state
      codC.rOc().rPerturbation().rLShell(0)(0,0) += 1.0/4.0;
      codC.rOc().rPerturbation().rLShell(0)(1,0) += -1.0/2.0;
   }

   template <typename TSimTraits> void SakurabaState<TSimTraits>::setSpecMagnetic(typename SakurabaState<TSimTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      if(PERTURB_MAGNETIC)
      {
         for(int l=1; l < pTrunc->local()->spec()->nL()/2; ++l)
         {
            // Set some perturbation random energy in Toroidal component
            magB.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l)/2,1)).setRandom();
            magB.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(RandomState<TSimTraits>::PERTURBATION_AMPLITUDE,0.0);

            // Set some perturbation random energy in Poloidal component
            magB.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l)/2,1)).setRandom();
            magB.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(RandomState<TSimTraits>::PERTURBATION_AMPLITUDE,0.0);
         }
      } else
      {
         for(int l=1; l < pTrunc->local()->spec()->nL()/2; ++l)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).setConstant(EPMComplex(0.0,0.0));
            magB.rOc().rPerturbation().rPol().rLShell(l).setConstant(EPMComplex(0.0,0.0));
         }
      }
   }

   template <typename TSimTraits> void SakurabaState<TSimTraits>::setSpecVelocity(typename SakurabaState<TSimTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/2; ++l)
      {
         velV.rOc().rPerturbation().rTor().rLShell(l).setConstant(EPMComplex(0.0,0.0));
         velV.rOc().rPerturbation().rPol().rLShell(l).setConstant(EPMComplex(0.0,0.0));
      }
   }

   template <typename TSimTraits> SakurabaState<TSimTraits>::SakurabaState()
   {
   }

}

#endif // SAKURABASTATE_HPP
