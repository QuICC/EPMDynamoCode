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

   template <typename TGenTraits> const EPMFloat BenchmarkState<TGenTraits>::PERTURBATION_AMPLITUDE = 0.0;

   template <typename TGenTraits> const int BenchmarkState<TGenTraits>::PERTURBATION_LRATIO = 1;

   template <typename TGenTraits> const int BenchmarkState<TGenTraits>::PERTURBATION_NRATIO = 2;

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setRTPCodensity(typename BenchmarkState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      Array sPh = pTrunc->sim()->hoz()->phGrid().array().sin();
      Array cPh = pTrunc->sim()->hoz()->phGrid().array().cos();
      Array s2Ph = (2.*pTrunc->sim()->hoz()->phGrid()).array().sin();
      Array c2Ph = (2.*pTrunc->sim()->hoz()->phGrid()).array().cos();
      Array s3Ph = (3.*pTrunc->sim()->hoz()->phGrid()).array().sin();
      Array c3Ph = (3.*pTrunc->sim()->hoz()->phGrid()).array().cos();
      Array s4Ph = (4.*pTrunc->sim()->hoz()->phGrid()).array().sin();
      Array c4Ph = (4.*pTrunc->sim()->hoz()->phGrid()).array().cos();

      EPMFloat r;
      EPMFloat r2;
      EPMFloat r3;
      EPMFloat r4;

      EPMFloat ampl = 1e-8;
      EPMFloat norm = (1./8.)*std::sqrt(35./MathConstants::PI);

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         r = pTrunc->local()->rtp()->radGrid(n);
         r2 = r*r;
         r3 = r2*r;
         r4 = r2*r2;

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            for(int ph=0; ph < sPh.size(); ++ph)
            {
               /////////////////////////////////////////////////////////////////////////////////////////////////////
               // 
               // FSBenchmark: Thermal convection at E=3e-4, Ra=95, Pr=1, no-slip, fixed temperature
               //
               // Extact approximation
               //codC.rOc().rRTP().rShell(n).col(th)(ph) = 0.48397401786548866-0.4810054713494751*r2;
               
               // Conducting state
               codC.rOc().rRTP().rShell(n).col(th)(ph) = 0.5-0.5*r2;

               // Y_2^2 perturbation
               //codC.rOc().rRTP().rShell(n).col(th)(ph) += ampl*norm*r2*(1.0 - r2)*s2Ph(ph)*std::pow(pTrunc->local()->rtp()->sTh(th,n),2);
               //codC.rOc().rRTP().rShell(n).col(th)(ph) += ampl*norm*r2*(1.0 - r2)*c2Ph(ph)*std::pow(pTrunc->local()->rtp()->sTh(th,n),2);
               
               // Y_3^3 perturbation
               //codC.rOc().rRTP().rShell(n).col(th)(ph) += ampl*norm*r3*(1.0 - r2)*s3Ph(ph)*std::pow(pTrunc->local()->rtp()->sTh(th,n),3);
               //codC.rOc().rRTP().rShell(n).col(th)(ph) += ampl*norm*r3*(1.0 - r2)*c3Ph(ph)*std::pow(pTrunc->local()->rtp()->sTh(th,n),3);
               
               // Y_4^4 perturbation
               codC.rOc().rRTP().rShell(n).col(th)(ph) += ampl*norm*r4*(1.0 - r2)*s4Ph(ph)*std::pow(pTrunc->local()->rtp()->sTh(th,n),4);
               codC.rOc().rRTP().rShell(n).col(th)(ph) += ampl*norm*r4*(1.0 - r2)*c4Ph(ph)*std::pow(pTrunc->local()->rtp()->sTh(th,n),4);

               //
               /////////////////////////////////////////////////////////////////////////////////////////////////////
               
               /////////////////////////////////////////////////////////////////////////////////////////////////////
               // 
               // FSBenchmark: Thermal convection at E=3e-4, Ra=95, Pr=1, stress-free, fixed temperature
               //
               // Extact approximation
               //codC.rOc().rRTP().rShell(n).col(th)(ph) = 0.49386130276515167-0.4925442767121938*r2;

               // Coarse approximation
               //codC.rOc().rRTP().rShell(n).col(th)(ph) = 0.5-0.5*r2;
               //
               /////////////////////////////////////////////////////////////////////////////////////////////////////
            }
         }
      }
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setRTPMagnetic(typename BenchmarkState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      Array sPh = pTrunc->sim()->hoz()->phGrid().array().sin();
      Array cPh = pTrunc->sim()->hoz()->phGrid().array().cos();

      EPMFloat r;
      EPMFloat r2;
      EPMFloat r4;
      EPMFloat r6;

      EPMFloat tmp;
      EPMFloat sTh;
      EPMFloat s2Th;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
         
         r = pTrunc->local()->rtp()->radGrid(n);
         r2 = r*r;
         r4 = r2*r2;
         r6 = r4*r2;

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            tmp = std::acos(pTrunc->local()->rtp()->cTh(th,n));
            sTh = pTrunc->local()->rtp()->sTh(th,n);
            s2Th = std::sin(2.0*tmp);

            for(int ph=0; ph < sPh.size(); ++ph)
            {
               // Radial component
               magB.rOc().rRTP().rR().rShell(n).col(th)(ph) = 0.0;

               // Theta component
               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = ((-3./2.)*r)*(-1. + 4.*r2 - 6.*r4 + 3.*r6)*(cPh(ph) + sPh(ph));

               // Phi component
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (3.*r)*(2. - 5.*r2 + 4.*r4)*sTh;
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += 2.*(1. - 3.*r2 + 3.*r4)*(cPh(ph) - sPh(ph));
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) *= (-3./4.)*r*(-1. + r2)*pTrunc->local()->rtp()->cTh(th,n);
            }
         }
      }
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setRTPVelocity(typename BenchmarkState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      EPMFloat r;
      EPMFloat r2;
      EPMFloat r4;
      EPMFloat r6;
      EPMFloat r8;

      EPMFloat sTh;
      EPMFloat cTh;
      EPMFloat s2Th;
      EPMFloat c2Th;
      EPMFloat s3Th;
      EPMFloat c3Th;
      EPMFloat s4Th;
      EPMFloat c4Th;

      EPMFloat sPh;
      EPMFloat cPh;

      EPMFloat tmp;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
         
         r = pTrunc->local()->rtp()->radGrid(n);
         r2 = r*r;
         r4 = r2*r2;
         r6 = r4*r2;
         r8 = r4*r4;

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            // Radial component
            velV.rOc().rRTP().rR().rShell(n).col(th).setConstant(0.0);

            sTh = pTrunc->local()->rtp()->sTh(th,n);
            cTh = pTrunc->local()->rtp()->cTh(th,n);
            tmp = std::acos(pTrunc->local()->rtp()->cTh(th,n));
            s2Th = std::sin(2.0*tmp);
            c2Th = std::cos(2.0*tmp);
            s3Th = std::sin(3.0*tmp);
            c3Th = std::cos(3.0*tmp);

            for(int ph=0; ph < pTrunc->sim()->hoz()->phGrid().size(); ++ph)
            {
                sPh = pTrunc->sim()->hoz()->phGrid().array().sin()(ph);
                cPh = pTrunc->sim()->hoz()->phGrid().array().cos()(ph);

                velV.rOc().rRTP().rR().rShell(n).col(th)(ph) = 0.0;

                velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = 3.*(-147. + 343.*r2 -217.*r4 + 29.*r6)*cPh;
                velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += 14.*(-9. - 125.*r2 + 39.*r4 + 27*r6)*sPh;
                velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) *= (-10./(7.*std::sqrt(3.)))*r2*cTh;

                velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = 7.*((43700. - 58113.*r2 - 15345.*r4 + 1881.*r6 + 20790.*r8)*sTh + 1485.*r2*(-9. + 115.*r2 - 167.*r4 + 70.*r6)*s3Th);
                velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += 528.*std::sqrt(3.)*c2Th*r*(14.*(-9 - 125.*r2 + 39.*r4 + 27.*r6)*cPh + 3.*(147. - 343.*r2 + 217.*r4 - 29.*r6)*sPh);
                velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) *= (-5./5544.)*r;
            }
         }
      }
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setSpecCodensity(typename BenchmarkState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      // Set some perturbation random energy
      for(int l=0; l < pTrunc->local()->spec()->nL()/BenchmarkState<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
//         codC.rOc().rPerturbation().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/BenchmarkState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
//         codC.rOc().rPerturbation().rLShell(l) *= EPMComplex(BenchmarkState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            codC.rOc().rPerturbation().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }

//      // Create basic state to generate basic -b r gradient (including normalisation factors)
//      codC.rOc().rPerturbation().rLShell(0)(0,0) += 0.25*(std::sqrt(MathConstants::PI/2.0));
//      codC.rOc().rPerturbation().rLShell(0)(1,0) += -0.5*(std::sqrt(MathConstants::PI)/4.0);
//
//      // Include perturbation of degree and order 4
//      EPMFloat amp = 0.0;
//      codC.rOc().rPerturbation().rLShell(4)(0,4) += 0.25*(std::sqrt(35.0*MathConstants::PI)/64.0)*amp;
//      codC.rOc().rPerturbation().rLShell(4)(1,4) += -0.5*(std::sqrt(105.0*MathConstants::PI/2.0)/64.0)*amp;

      // Use approximate solution (> 1e-1)
      codC.rOc().rPerturbation().rLShell(0)(0,0) += 0.2748;
      codC.rOc().rPerturbation().rLShell(0)(1,0) += -0.1834;
   }

   template <typename TGenTraits> void BenchmarkState<TGenTraits>::setSpecMagnetic(typename BenchmarkState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL(); ++l)
      {
         // Make sure toroidal components are initialised to zero
         magB.rOc().rPerturbation().rTor().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Make sure poloidal components are initialised to zero
         magB.rOc().rPerturbation().rPol().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }

      EPMFloat ampT = 3.0;
      EPMFloat ampP = 0.0;

      magB.rOc().rPerturbation().rTor().rLShell(1)(0,1) = ampT*EPMComplex(0.0242328,0.0242328);
      magB.rOc().rPerturbation().rTor().rLShell(1)(1,1) = ampT*EPMComplex(-0.0242328,-0.0242328);
      magB.rOc().rPerturbation().rTor().rLShell(1)(2,1) = ampT*EPMComplex(0.0103855,0.0103855);
      magB.rOc().rPerturbation().rTor().rLShell(1)(3,1) = ampT*EPMComplex(-0.0103855,-0.0103855);
      magB.rOc().rPerturbation().rTor().rLShell(2)(0,0) = ampT*EPMComplex(0.0199868,0.0);
      magB.rOc().rPerturbation().rTor().rLShell(2)(1,0) = ampT*EPMComplex(-0.0250275,0.0);
      magB.rOc().rPerturbation().rTor().rLShell(2)(2,0) = ampT*EPMComplex(0.0105328,0.0);
      magB.rOc().rPerturbation().rTor().rLShell(2)(3,0) = ampT*EPMComplex(-0.00785068,0.0);
      //magB.rOc().rPerturbation().rTor().rLShell(2)(0,1) = ampT*EPMComplex(-0.00239842,0.0000799474);
      //magB.rOc().rPerturbation().rTor().rLShell(2)(1,1) = ampT*EPMComplex(0.00321782,0.000393289);
      //magB.rOc().rPerturbation().rTor().rLShell(2)(2,1) = ampT*EPMComplex(0.000234062,-0.000257468);
      //magB.rOc().rPerturbation().rTor().rLShell(2)(3,1) = ampT*EPMComplex(-0.000785068,-0.00023552);

      magB.rOc().rPerturbation().rPol().rLShell(1)(0,0) = ampP*EPMComplex(0.00249251,0.0);
      magB.rOc().rPerturbation().rPol().rLShell(1)(1,0) = ampP*EPMComplex(-0.00249251,0.0);
      magB.rOc().rPerturbation().rPol().rLShell(1)(2,0) = ampP*EPMComplex(0.000830838,0.0);
      magB.rOc().rPerturbation().rPol().rLShell(1)(3,0) = ampP*EPMComplex(-0.000830838,0.0);
      magB.rOc().rPerturbation().rPol().rLShell(2)(0,1) = ampP*EPMComplex(0.000959369,-0.000159895);
      magB.rOc().rPerturbation().rPol().rLShell(2)(1,1) = ampP*EPMComplex(-0.000858086,0.000786578);
      magB.rOc().rPerturbation().rPol().rLShell(2)(2,1) = ampP*EPMComplex(-0.000842623,-0.00098306);
      magB.rOc().rPerturbation().rPol().rLShell(2)(3,1) = ampP*EPMComplex(0.000628054,0.000366365);
      magB.rOc().rPerturbation().rPol().rLShell(3)(0,0) = ampP*EPMComplex(0.00109472,0.0);
      magB.rOc().rPerturbation().rPol().rLShell(3)(1,0) = ampP*EPMComplex(-0.00277563,0.0);
      magB.rOc().rPerturbation().rPol().rLShell(3)(2,0) = ampP*EPMComplex(0.000201114,0.0);
      magB.rOc().rPerturbation().rPol().rLShell(3)(3,0) = ampP*EPMComplex(-0.000485185,0.0);

//      // Toroidal perturbation of degree 2
//      EPMFloat ampTor = 1.0;
//      int torL = 2;
//      EPMFloat torN = (1.0)/(1.0 + static_cast<EPMFloat>(torL));
//      EPMFloat torN0 = 0.5*torN;
//      EPMFloat torN1 = torN;
//      magB.rOc().rPerturbation().rTor().rLShell(torL)(0,0) += torN0*(std::sqrt(3.0*MathConstants::PI)/4.0)*ampTor;
//      magB.rOc().rPerturbation().rTor().rLShell(torL)(1,0) += -torN1*(std::sqrt(15.0*MathConstants::PI)/16.0)*ampTor;
//      magB.rOc().rPerturbation().rTor().rLShell(torL)(0,1) += torN0*(std::sqrt(3.0*MathConstants::PI)/4.0)*ampTor;
//      magB.rOc().rPerturbation().rTor().rLShell(torL)(1,1) += -torN1*(std::sqrt(15.0*MathConstants::PI)/16.0)*ampTor;
//      magB.rOc().rPerturbation().rTor().rLShell(torL)(0,2) += torN0*(std::sqrt(3.0*MathConstants::PI)/4.0)*ampTor;
//      magB.rOc().rPerturbation().rTor().rLShell(torL)(1,2) += -torN1*(std::sqrt(15.0*MathConstants::PI)/16.0)*ampTor;
//
//      // Poloidal perturbation of degree 1
//      EPMFloat ampPol = 1.0;
//      int polL = 1;
//      EPMFloat polN0 = (5.0 + 6.0*static_cast<EPMFloat>(polL))/(6.0 + 10.0*static_cast<EPMFloat>(polL) + 4.0*static_cast<EPMFloat>(polL)*static_cast<EPMFloat>(polL));
//      EPMFloat polN1 = (1.0 + 2.0*static_cast<EPMFloat>(polL))/(3.0 + 5.0*static_cast<EPMFloat>(polL) + 2.0*static_cast<EPMFloat>(polL)*static_cast<EPMFloat>(polL));
//      magB.rOc().rPerturbation().rPol().rLShell(polL)(0,0) += polN0*(std::sqrt(MathConstants::PI)/2.0)*ampPol;
//      magB.rOc().rPerturbation().rPol().rLShell(polL)(1,0) += -polN1*(std::sqrt(MathConstants::PI)/4.0)*ampPol;
////      magB.rOc().rPerturbation().rPol().rLShell(polL)(0,1) += polN0*(std::sqrt(MathConstants::PI)/2.0)*ampPol;
////      magB.rOc().rPerturbation().rPol().rLShell(polL)(1,1) += -polN1*(std::sqrt(MathConstants::PI)/4.0)*ampPol;
      
      // Use approximate toroidal solution (>1e-2)
      // l=1, m=1
//      magB.rOc().rPerturbation().rTor().rLShell(1)(2,1) += -0.0213*MathConstants::cI;
//      magB.rOc().rPerturbation().rTor().rLShell(1)(3,1) += 0.0173*MathConstants::cI;
//      // l=2, m=0
//      magB.rOc().rPerturbation().rTor().rLShell(2)(0,0) += 0.0429;
//      magB.rOc().rPerturbation().rTor().rLShell(2)(1,0) += -0.0445;
//      magB.rOc().rPerturbation().rTor().rLShell(2)(3,0) += -0.0201;
//      magB.rOc().rPerturbation().rTor().rLShell(2)(4,0) += 0.0162;
//      // l=3, m=1
//      magB.rOc().rPerturbation().rTor().rLShell(3)(1,1) += -0.0105;
//      // l=4, m=0
//      magB.rOc().rPerturbation().rTor().rLShell(4)(0,0) += -0.0126;
//      // l=5, m=1
//      magB.rOc().rPerturbation().rTor().rLShell(5)(1,1) += -0.0116*MathConstants::cI;

      // Use approximate poloidal solution (>1e-2)
      // l=1, m=0
//      magB.rOc().rPerturbation().rPol().rLShell(1)(0,0) += 0.0190;
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

      EPMFloat ampT = 1.0;
      EPMFloat ampP = 0.0;

      velV.rOc().rPerturbation().rTor().rLShell(1)(0,0) = ampT*EPMComplex(3.71709,0.0);
      velV.rOc().rPerturbation().rTor().rLShell(1)(1,0) = ampT*EPMComplex(69.5827,0.0);
      velV.rOc().rPerturbation().rTor().rLShell(1)(2,0) = ampT*EPMComplex(-14.1935,0.0);
      velV.rOc().rPerturbation().rTor().rLShell(1)(3,0) = ampT*EPMComplex(-4.50037,0.0);
      velV.rOc().rPerturbation().rTor().rLShell(2)(0,1) = ampT*EPMComplex(171.434,0.199868);
      velV.rOc().rPerturbation().rTor().rLShell(2)(1,1) = ampT*EPMComplex(13.2884,-8.67024);
      velV.rOc().rPerturbation().rTor().rLShell(2)(2,1) = ampT*EPMComplex(-9.94764,3.77007);
      velV.rOc().rPerturbation().rTor().rLShell(2)(3,1) = ampT*EPMComplex(-0.706561,-0.162621);
      velV.rOc().rPerturbation().rTor().rLShell(3)(0,0) = ampT*EPMComplex(-37.2863,0.0);
      velV.rOc().rPerturbation().rTor().rLShell(3)(1,0) = ampT*EPMComplex(6.08789,0.0);
      velV.rOc().rPerturbation().rTor().rLShell(3)(2,0) = ampT*EPMComplex(2.42677,0.0);
      velV.rOc().rPerturbation().rTor().rLShell(3)(3,0) = ampT*EPMComplex(-1.69815,0.0);

      velV.rOc().rPerturbation().rPol().rLShell(1)(0,1) = ampP*EPMComplex(-0.247713,-0.721598);
      velV.rOc().rPerturbation().rPol().rLShell(1)(1,1) = ampP*EPMComplex(1.16471,0.580048);
      velV.rOc().rPerturbation().rPol().rLShell(1)(2,1) = ampP*EPMComplex(-1.18317,0.213864);
      velV.rOc().rPerturbation().rPol().rLShell(1)(3,1) = ampP*EPMComplex(0.266176,-0.0723137);
      velV.rOc().rPerturbation().rPol().rLShell(2)(0,0) = ampP*EPMComplex(-0.0470801,0.0);
      velV.rOc().rPerturbation().rPol().rLShell(2)(1,0) = ampP*EPMComplex(0.0297946,0.0);
      velV.rOc().rPerturbation().rPol().rLShell(2)(2,0) = ampP*EPMComplex(0.0351093,0.0);
      velV.rOc().rPerturbation().rPol().rLShell(2)(3,0) = ampP*EPMComplex(-0.0119214,0.0);
      velV.rOc().rPerturbation().rPol().rLShell(2)(0,2) = ampP*EPMComplex(-1.53745,-0.706202);
      velV.rOc().rPerturbation().rPol().rLShell(2)(1,2) = ampP*EPMComplex(1.9747,1.12624);
      velV.rOc().rPerturbation().rPol().rLShell(2)(2,2) = ampP*EPMComplex(-0.234062,-0.409609);
      velV.rOc().rPerturbation().rPol().rLShell(2)(3,2) = ampP*EPMComplex(-0.0265715,0.0671669);
      velV.rOc().rPerturbation().rPol().rLShell(3)(0,1) = ampP*EPMComplex(-0.12168,0.0207003);
      velV.rOc().rPerturbation().rPol().rLShell(3)(1,1) = ampP*EPMComplex(0.203546,0.0);
      velV.rOc().rPerturbation().rPol().rLShell(3)(2,1) = ampP*EPMComplex(-0.072401,-0.0375412);
      velV.rOc().rPerturbation().rPol().rLShell(3)(3,1) = ampP*EPMComplex(0.0108799,0.0123502);

//      // Use approximate poloidal solution (>1e0)
//      // l=1, m=1
//      velV.rOc().rPerturbation().rPol().rLShell(1)(0,1) += EPMComplex(1.1045,-0.1641);
//      velV.rOc().rPerturbation().rPol().rLShell(1)(1,1) += EPMComplex(-1.8535,-0.4113);
//      velV.rOc().rPerturbation().rPol().rLShell(1)(2,1) += EPMComplex(0.7476,1.0009);
//      velV.rOc().rPerturbation().rPol().rLShell(1)(3,1) += EPMComplex(0.0,-0.5084);
//      // l=2, m=0
//      velV.rOc().rPerturbation().rPol().rLShell(2)(0,0) += EPMComplex(-0.1064,0.0);
//      velV.rOc().rPerturbation().rPol().rLShell(2)(1,0) += EPMComplex(0.1126,0.0);
//      // l=2, m=2
//      velV.rOc().rPerturbation().rPol().rLShell(2)(0,2) += EPMComplex(-1.6291,0.5112);
//      velV.rOc().rPerturbation().rPol().rLShell(2)(1,2) += EPMComplex(1.7500,0.0);
//      velV.rOc().rPerturbation().rPol().rLShell(2)(2,2) += EPMComplex(0.3940,-0.6714);
//      velV.rOc().rPerturbation().rPol().rLShell(2)(3,2) += EPMComplex(-0.4246,0.0);
//      velV.rOc().rPerturbation().rPol().rLShell(2)(4,2) += EPMComplex(0.1256,0.0);
//      // l=3, m=1
//      velV.rOc().rPerturbation().rPol().rLShell(3)(1,1) += EPMComplex(0.0,-0.1559);
//      velV.rOc().rPerturbation().rPol().rLShell(3)(2,1) += EPMComplex(0.0,0.1319);
//      // l=3, m=3
//      velV.rOc().rPerturbation().rPol().rLShell(3)(0,3) += EPMComplex(0.7419,0.0);
//      velV.rOc().rPerturbation().rPol().rLShell(3)(1,3) += EPMComplex(-0.7475,-0.2183);
//      velV.rOc().rPerturbation().rPol().rLShell(3)(2,3) += EPMComplex(-0.3438,0.1713);
//      velV.rOc().rPerturbation().rPol().rLShell(3)(3,3) += EPMComplex(0.2255,0.2397);
//      velV.rOc().rPerturbation().rPol().rLShell(3)(5,3) += EPMComplex(0.0,-0.1319);
//      // l=4, m=2
//      velV.rOc().rPerturbation().rPol().rLShell(4)(1,2) += EPMComplex(-0.1021,0.1360);
//      velV.rOc().rPerturbation().rPol().rLShell(4)(2,2) += EPMComplex(0.1432,0.0);
//      // l=4, m=4
//      velV.rOc().rPerturbation().rPol().rLShell(4)(0,4) += EPMComplex(0.0,-0.3576);
//      velV.rOc().rPerturbation().rPol().rLShell(4)(1,4) += EPMComplex(-0.2379,0.3798);
//      velV.rOc().rPerturbation().rPol().rLShell(4)(2,4) += EPMComplex(0.1888,0.2094);
//      velV.rOc().rPerturbation().rPol().rLShell(4)(3,4) += EPMComplex(0.0,-0.1886);
//      velV.rOc().rPerturbation().rPol().rLShell(4)(4,4) += EPMComplex(-0.1134,0.0);
//      // l=5, m=5
//      velV.rOc().rPerturbation().rPol().rLShell(5)(1,5) += EPMComplex(0.0,0.1778);
//      velV.rOc().rPerturbation().rPol().rLShell(5)(2,5) += EPMComplex(0.0,-0.1302);
//      velV.rOc().rPerturbation().rPol().rLShell(5)(3,5) += EPMComplex(-0.1129,0.0);
   }

   template <typename TGenTraits> BenchmarkState<TGenTraits>::BenchmarkState()
   {
   }
}

#endif // BENCHMARKSTATE_HPP
