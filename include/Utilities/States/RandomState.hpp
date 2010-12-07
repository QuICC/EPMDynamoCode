/** \file RandomState.hpp
 *  \brief Initial state definition function for a random perturbation state
 */

#ifndef RANDOMSTATE_HPP
#define RANDOMSTATE_HPP

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
   class RandomTraits
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
    * \brief Initial state definition function for a random perturbation state
    *
    * \epmBug State generation is not aware of parallelisation
    */
   template <typename TGenTraits> class RandomState
   {
      public:
         /// Typdef for the StateTraits type
         typedef RandomTraits  StateTraits;

         /**
          * @brief Perturbation amplitude
          */
         static const EPMFloat PERTURBATION_AMPLITUDE;

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
         RandomState();

         /**
          * @brief Destructor
          */
         virtual ~RandomState() {};
   };

   template <typename TGenTraits> const EPMFloat RandomState<TGenTraits>::PERTURBATION_AMPLITUDE = 1.0e-5;

   template <typename TGenTraits> void RandomState<TGenTraits>::setRTPCodensity(typename RandomState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void RandomState<TGenTraits>::setRTPMagnetic(typename RandomState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void RandomState<TGenTraits>::setRTPVelocity(typename RandomState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
      }
   }

   template <typename TGenTraits> void RandomState<TGenTraits>::setSpecCodensity(typename RandomState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      // Set some perturbation random energy
      for(int l=0; l < pTrunc->local()->spec()->nL()/2; ++l)
      {
         codC.rOc().rPerturbation().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         codC.rOc().rPerturbation().rLShell(l) *= EPMComplex(RandomState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            codC.rOc().rPerturbation().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }

      // Create basic state
      codC.rOc().rPerturbation().rLShell(0)(0,0) += 1.0/4.0;
      codC.rOc().rPerturbation().rLShell(0)(1,0) += -1.0/2.0;
   }

   template <typename TGenTraits> void RandomState<TGenTraits>::setSpecMagnetic(typename RandomState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/2; ++l)
      {
         // Set some perturbation random energy in Toroidal component
         magB.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(RandomState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Set some perturbation random energy in Poloidal component
         magB.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(RandomState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }
   }

   template <typename TGenTraits> void RandomState<TGenTraits>::setSpecVelocity(typename RandomState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/2; ++l)
      {
         // Set some perturbation random energy in Toroidal component
         velV.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         velV.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(RandomState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Set some perturbation random energy in Poloidal component
         velV.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         velV.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(RandomState<TGenTraits>::PERTURBATION_AMPLITUDE, 0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            velV.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            velV.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }
   }

   template <typename TGenTraits> RandomState<TGenTraits>::RandomState()
   {
   }
}

#endif // RANDOMSTATE_HPP
