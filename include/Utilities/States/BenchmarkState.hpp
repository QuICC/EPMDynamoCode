/** \file BenchmarkState.hpp
 *  \brief Initial state definition(s) for the dynamo benchmark(s)
 */

#ifndef BENCHMARKSTATE_HPP
#define BENCHMARKSTATE_HPP

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
    * @brief Traits of benchmark intial state generator
    */
   class BenchmarkTraits
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
    * \brief Initial state definition(s) for the dynamo benchmark(s)
    *
    * \epmBug State generation is not aware of parallelisation
    */
   template <typename TGenTraits> class BenchmarkState
   {
      public:
         /// Typdef for the StateTraits type
         typedef BenchmarkTraits  StateTraits;

         /**
          * @brief Benchmark case state
          */
         static const int BENCHMARK_CASE;

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
         BenchmarkState();

         /**
          * @brief Destructor
          */
         virtual ~BenchmarkState() {};
   };

   template <typename TGenTraits> const int BenchmarkState<TGenTraits>::BENCHMARK_CASE = 0;

   template <typename TGenTraits> const EPMFloat BenchmarkState<TGenTraits>::PERTURBATION_AMPLITUDE = 1.0e-10;

   template <typename TGenTraits> const int BenchmarkState<TGenTraits>::PERTURBATION_LRATIO = 1;

   template <typename TGenTraits> const int BenchmarkState<TGenTraits>::PERTURBATION_NRATIO = 2;

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setRTPCodensity(typename BenchmarkState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setRTPMagnetic(typename BenchmarkState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      Array sPh = pTrunc->sim()->hoz()->phGrid().array().sin();
      Array cPh = pTrunc->sim()->hoz()->phGrid().array().cos();

      EPMFloat radius;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
         
         radius = pTrunc->local()->rtp()->radGrid(n);

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            // Radial component
            magB.rOc().rRTP().rR().rShell(n).col(th).setConstant(1.25*(4.0-3.0*radius)*pTrunc->local()->rtp()->cTh(th, n));

            // Theta component
            magB.rOc().rRTP().rTheta().rShell(n).col(th).setConstant(0.6125*(9.0*radius-8.0)*pTrunc->local()->rtp()->sTh(th, n));

            // Phi component
            magB.rOc().rRTP().rPhi().rShell(n).col(th).setConstant(5.0*std::sin(std::acos(-1.0)*radius)*2.0*pTrunc->local()->rtp()->sTh(th, n)*pTrunc->local()->rtp()->cTh(th, n));
         }
      }
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setRTPVelocity(typename BenchmarkState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
      }
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setSpecCodensity(typename BenchmarkState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      // Set some perturbation random energy
      for(int l=0; l < pTrunc->local()->spec()->nL()/BenchmarkState<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         codC.rOc().rPerturbation().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/BenchmarkState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         codC.rOc().rPerturbation().rLShell(l) *= EPMComplex(BenchmarkState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            codC.rOc().rPerturbation().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }

      // Create basic state to generate basic -b r gradient (including normalisation factors)
      codC.rOc().rPerturbation().rLShell(0)(0,0) += 0.25*(std::sqrt(MathConstants::PI/2.0));
      codC.rOc().rPerturbation().rLShell(0)(1,0) += -0.5*(std::sqrt(MathConstants::PI)/4.0);

      // Include perturbation of degree and order 4
      EPMFloat amp = 0.0;
      codC.rOc().rPerturbation().rLShell(4)(0,4) += 0.25*(std::sqrt(35.0*MathConstants::PI)/64.0)*amp;
      codC.rOc().rPerturbation().rLShell(4)(1,4) += -0.5*(std::sqrt(105.0*MathConstants::PI/2.0)/64.0)*amp;
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setSpecMagnetic(typename BenchmarkState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL(); ++l)
      {
         // Set some perturbation random energy in Toroidal component
         magB.rOc().rPerturbation().rTor().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Set some perturbation random energy in Poloidal component
         magB.rOc().rPerturbation().rPol().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }

      // Toroidal perturbation of degree 2
      EPMFloat ampTor = 1.0;
      int torL = 2;
      EPMFloat torN = (1.0)/(1.0 + static_cast<EPMFloat>(torL));
      EPMFloat torN0 = 0.5*torN;
      EPMFloat torN1 = torN;
      magB.rOc().rPerturbation().rTor().rLShell(torL)(0,0) += torN0*(std::sqrt(3.0*MathConstants::PI)/4.0)*ampTor;
      magB.rOc().rPerturbation().rTor().rLShell(torL)(1,0) += -torN1*(std::sqrt(15.0*MathConstants::PI)/16.0)*ampTor;
      magB.rOc().rPerturbation().rTor().rLShell(torL)(0,1) += torN0*(std::sqrt(3.0*MathConstants::PI)/4.0)*ampTor;
      magB.rOc().rPerturbation().rTor().rLShell(torL)(1,1) += -torN1*(std::sqrt(15.0*MathConstants::PI)/16.0)*ampTor;
      magB.rOc().rPerturbation().rTor().rLShell(torL)(0,2) += torN0*(std::sqrt(3.0*MathConstants::PI)/4.0)*ampTor;
      magB.rOc().rPerturbation().rTor().rLShell(torL)(1,2) += -torN1*(std::sqrt(15.0*MathConstants::PI)/16.0)*ampTor;

      // Poloidal perturbation of degree 1
      EPMFloat ampPol = 1.0;
      int polL = 1;
      EPMFloat polN0 = (5.0 + 6.0*static_cast<EPMFloat>(polL))/(6.0 + 10.0*static_cast<EPMFloat>(polL) + 4.0*static_cast<EPMFloat>(polL)*static_cast<EPMFloat>(polL));
      EPMFloat polN1 = (1.0 + 2.0*static_cast<EPMFloat>(polL))/(3.0 + 5.0*static_cast<EPMFloat>(polL) + 2.0*static_cast<EPMFloat>(polL)*static_cast<EPMFloat>(polL));
      magB.rOc().rPerturbation().rPol().rLShell(polL)(0,0) += polN0*(std::sqrt(MathConstants::PI)/2.0)*ampPol;
      magB.rOc().rPerturbation().rPol().rLShell(polL)(1,0) += -polN1*(std::sqrt(MathConstants::PI)/4.0)*ampPol;
//      magB.rOc().rPerturbation().rPol().rLShell(polL)(0,1) += polN0*(std::sqrt(MathConstants::PI)/2.0)*ampPol;
//      magB.rOc().rPerturbation().rPol().rLShell(polL)(1,1) += -polN1*(std::sqrt(MathConstants::PI)/4.0)*ampPol;
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setSpecVelocity(typename BenchmarkState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL(); ++l)
      {
         // Set some perturbation random energy in Toroidal component
         velV.rOc().rPerturbation().rTor().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Set some perturbation random energy in Poloidal component
         velV.rOc().rPerturbation().rPol().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            velV.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            velV.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }
   }

   template <typename TGenTraits> BenchmarkState<TGenTraits>::BenchmarkState()
   {
   }
}

#endif // BENCHMARKSTATE_HPP
