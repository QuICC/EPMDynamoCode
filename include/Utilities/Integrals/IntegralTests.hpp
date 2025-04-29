/** \file IntegralTests.hpp
 *  \brief Initial state definition function for a random perturbation state
 */

#ifndef INTEGRALTESTS_HPP
#define INTEGRALTESTS_HPP

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
   class TestsTraits
   {
      public:
         /// Requires RTP Codensity computations
         static const bool UseRTPCodensity = true;

         /// Requires spectral Codensity computations
         static const bool UseSpecCodensity = false;

         /// Requires Codensity gradient computations
         static const bool UseCodensityGrad = false;

         /// Requires RTP Magnetic computations
         static const bool UseRTPMagnetic = true;

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
    */
   template <typename TGenTraits> class IntegralTests
   {
      public:
         /// Typdef for the IntegralTraits type
         typedef TestsTraits  IntegralTraits;

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
         IntegralTests();

         /**
          * @brief Destructor
          */
         virtual ~IntegralTests() {};
   };

   template <typename TGenTraits> const EPMFloat IntegralTests<TGenTraits>::PERTURBATION_AMPLITUDE = 1.0e-10;

   template <typename TGenTraits> const int IntegralTests<TGenTraits>::PERTURBATION_LRATIO = 2;

   template <typename TGenTraits> const int IntegralTests<TGenTraits>::PERTURBATION_NRATIO = 4;

   template <typename TGenTraits> void IntegralTests<TGenTraits>::setRTPCodensity(typename IntegralTests<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      int rPow = 1;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(std::pow(pTrunc->sim()->rad()->radGrid()(n), static_cast<EPMFloat>(rPow)/2.0));
      }
   }

   template <typename TGenTraits> void IntegralTests<TGenTraits>::setRTPMagnetic(typename IntegralTests<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      Array sPh = pTrunc->sim()->hoz()->phGrid().array().sin();
      Array cPh = pTrunc->sim()->hoz()->phGrid().array().cos();

      EPMFloat r;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         r = pTrunc->local()->rtp()->radGrid(n);

         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            magB.rOc().rRTP().rR().rShell(n).col(th).setConstant(2.0*pTrunc->local()->rtp()->cTh(th, n));
            magB.rOc().rRTP().rTheta().rShell(n).col(th).setConstant(-2.0*pTrunc->local()->rtp()->sTh(th, n));
            magB.rOc().rRTP().rPhi().rShell(n).col(th).setConstant(r*pTrunc->local()->rtp()->sTh(th, n));
         }
      }
   }

   template <typename TGenTraits> void IntegralTests<TGenTraits>::setRTPVelocity(typename IntegralTests<TGenTraits>::Velocity &velV)
   {

      SmartTruncation pTrunc = velV.oc().trunc();

      EPMFloat r;
      EPMFloat r2;
      EPMFloat r4;
      EPMFloat r8;

      EPMFloat sTh;
      EPMFloat cTh;

      EPMFloat sPh;
      EPMFloat cPh;
      EPMFloat s2Ph;
      EPMFloat c2Ph;
      EPMFloat s3Ph;
      EPMFloat c3Ph;

      EPMFloat tmp;

      // Initialise field with   Tor = r^3(1+r^2+r^4+r^6)Y_3^3
      // ... and                 Pol = r^3(1+r^2+r^4+r^6)Y_3^3

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);

         r = pTrunc->local()->rtp()->radGrid(n);
         r2 = r*r;
         r4 = r2*r2;
         r8 = r4*r4;

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            sTh = pTrunc->local()->rtp()->sTh(th,n);
            cTh = pTrunc->local()->rtp()->cTh(th,n);

            for(int ph=0; ph < pTrunc->sim()->hoz()->phGrid().size(); ++ph)
            {
               tmp = pTrunc->sim()->hoz()->phGrid()(ph);
               sPh = std::sin(tmp);
               cPh = std::cos(tmp);
               s2Ph = std::sin(2.0*tmp);
               c2Ph = std::cos(2.0*tmp);
               s3Ph = std::sin(3.0*tmp);
               c3Ph = std::cos(3.0*tmp);

               velV.rOc().rRTP().rR().rShell(n).col(th)(ph) = -3.0*std::sqrt(35/MathConstants::PI)*r2*(1.0+r2+r4+r8)*sTh*sTh*sTh*c3Ph;

               velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = (3.0/4.0)*std::sqrt(35/MathConstants::PI)*r2*sTh*sTh*((-2.0)*(2.0+3.0*r2+4.0*r4+6.0*r8)*cTh*c3Ph+r*(1.0+r2+r4+r8)*s3Ph);

               velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (3.0/4.0)*std::sqrt(35/MathConstants::PI)*r2*sTh*sTh*(r*(1.0+r2+r4+r8)*cTh*c3Ph+2.0*(2.0+3.0*r2+4.0*r4+6.0*r8)*s3Ph);
            }
         }
      }
   }

   template <typename TGenTraits> void IntegralTests<TGenTraits>::setSpecCodensity(typename IntegralTests<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      // Set some perturbation random energy
      for(int l=0; l < pTrunc->local()->spec()->nL()/IntegralTests<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         codC.rOc().rPerturbation().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/IntegralTests<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         codC.rOc().rPerturbation().rLShell(l) *= EPMComplex(IntegralTests<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            codC.rOc().rPerturbation().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }

      // Create basic state
      //codC.rOc().rPerturbation().rLShell(0)(0,0) += 1.0/4.0;
      codC.rOc().rPerturbation().rLShell(0)(0,0) += 0.313329;
      //codC.rOc().rPerturbation().rLShell(0)(1,0) += -1.0/2.0;
      codC.rOc().rPerturbation().rLShell(0)(1,0) += -0.221557;
   }

   template <typename TGenTraits> void IntegralTests<TGenTraits>::setSpecMagnetic(typename IntegralTests<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/IntegralTests<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         // Set some perturbation random energy in Toroidal component
         magB.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/IntegralTests<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(IntegralTests<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Set some perturbation random energy in Poloidal component
         magB.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/IntegralTests<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(IntegralTests<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }
   }

   template <typename TGenTraits> void IntegralTests<TGenTraits>::setSpecVelocity(typename IntegralTests<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/IntegralTests<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         // Set some perturbation random energy in Toroidal component
         velV.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/IntegralTests<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         velV.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(IntegralTests<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Set some perturbation random energy in Poloidal component
         velV.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/IntegralTests<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         velV.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(IntegralTests<TGenTraits>::PERTURBATION_AMPLITUDE, 0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            velV.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            velV.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }
   }

   template <typename TGenTraits> IntegralTests<TGenTraits>::IntegralTests()
   {
   }
}

#endif // INTEGRALTESTS_HPP
