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
      EPMFloat ampl = 1e-5;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
         
         r = pTrunc->local()->rtp()->radGrid(n);
         r2 = r*r;

         for(int th=0; th < pTrunc->local()->rtp()->nTh(n); ++th)
         {
            for(int ph=0; ph < sPh.size(); ++ph)
            {
               // Radial component
               magB.rOc().rRTP().rR().rShell(n).col(th)(ph) = ampl*std::sqrt(3./MathConstants::PI)*(1.-3./5.*r2)*pTrunc->local()->rtp()->cTh(th,n);
               magB.rOc().rRTP().rR().rShell(n).col(th)(ph) += ampl*1./5.*std::sqrt(3./(2.*MathConstants::PI))*(-5.+3.*r2)*pTrunc->local()->rtp()->sTh(th,n)*cPh(ph);
               magB.rOc().rRTP().rR().rShell(n).col(th)(ph) += ampl*1./5.*std::sqrt(3./(2.*MathConstants::PI))*(-5.+3.*r2)*pTrunc->local()->rtp()->sTh(th,n)*sPh(ph);

               // Theta component
               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = ampl*1./5.*std::sqrt(3./MathConstants::PI)*(-5.+6.*r2)*pTrunc->local()->rtp()->sTh(th,n);
               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ampl*1./5.*std::sqrt(3./(2.*MathConstants::PI))*(-5.+6.*r2)*pTrunc->local()->rtp()->cTh(th,n)*cPh(ph);
               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ampl*1./5.*std::sqrt(3./(2.*MathConstants::PI))*(-5.+6.*r2)*pTrunc->local()->rtp()->cTh(th,n)*sPh(ph);
               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += -ampl*1./2.*std::sqrt(3./(2.*MathConstants::PI))*r*(-1.+r2)*sPh(ph);
               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ampl*1./2.*std::sqrt(3./(2.*MathConstants::PI))*r*(-1.+r2)*cPh(ph);

               // Phi component
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = -ampl*1./2.*std::sqrt(3./MathConstants::PI)*r*(-1.+r2)*pTrunc->local()->rtp()->sTh(th,n);
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ampl*1./5.*std::sqrt(3./(2.*MathConstants::PI))*(-5.+6.*r2)*cPh(ph);
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ampl*1./5.*std::sqrt(3./(2.*MathConstants::PI))*(-5.+6.*r2)*sPh(ph);
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += -ampl*1./2.*std::sqrt(3./(2.*MathConstants::PI))*r*(-1.+r2)*pTrunc->local()->rtp()->cTh(th,n)*cPh(ph);
               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += -ampl*1./2.*std::sqrt(3./(2.*MathConstants::PI))*r*(-1.+r2)*pTrunc->local()->rtp()->cTh(th,n)*sPh(ph);

//               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = (-1.09974 + 4.39897*std::pow(radius,2.0) - 3.51917*std::pow(radius,4.0))*cPh(ph);
//               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += (-0.329923 + 0.4398978*std::pow(radius,2.0))*sPh(ph);
//               magB.rOc().rRTP().rTheta().rShell(n).col(th)(ph) *= radius;
//
//
//               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (-0.329923 + 0.4398978*std::pow(radius,2.0))*cPh(ph);
//               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += (1.09974 - 4.39897*std::pow(radius,2.0) + 3.51917*std::pow(radius,4.0))*sPh(ph);
//               magB.rOc().rRTP().rPhi().rShell(n).col(th)(ph) *= radius*pTrunc->local()->rtp()->cTh(th,n);
            }
         }
      }
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
      EPMFloat s2Ph;
      EPMFloat c2Ph;
      EPMFloat s3Ph;
      EPMFloat c3Ph;
      EPMFloat s4Ph;
      EPMFloat c4Ph;

      EPMFloat tmp;

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
//         velV.rOc().rRTP().rR().rShell(n).setConstant(1e-16);
//         velV.rOc().rRTP().rTheta().rShell(n).setConstant(1e-16);
//         velV.rOc().rRTP().rPhi().rShell(n).setConstant(1e-16);
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
            s4Th = std::sin(4.0*tmp);
            c4Th = std::cos(4.0*tmp);

            for(int ph=0; ph < pTrunc->sim()->hoz()->phGrid().size(); ++ph)
            {
                sPh = pTrunc->sim()->hoz()->phGrid().array().sin()(ph);
                cPh = pTrunc->sim()->hoz()->phGrid().array().cos()(ph);
                s2Ph = (2.0*pTrunc->sim()->hoz()->phGrid().array()).sin()(ph);
                c2Ph = (2.0*pTrunc->sim()->hoz()->phGrid().array()).cos()(ph);
                s3Ph = (3.0*pTrunc->sim()->hoz()->phGrid().array()).sin()(ph);
                c3Ph = (3.0*pTrunc->sim()->hoz()->phGrid().array()).cos()(ph);
                s4Ph = (4.0*pTrunc->sim()->hoz()->phGrid().array()).sin()(ph);
                c4Ph = (4.0*pTrunc->sim()->hoz()->phGrid().array()).cos()(ph);

                velV.rOc().rRTP().rR().rShell(n).col(th)(ph) += r*(-275.884 + 350.*r2 + (625./2.)*r4 - 325.*r6)*pTrunc->local()->rtp()->sTh(th,n);
                velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += r*(-275.884 + 350.*r2 + (625./2.)*r4 - 325.*r6)*pTrunc->local()->rtp()->sTh(th,n);
                velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += r*(-275.884 + 350.*r2 + (625./2.)*r4 - 325.*r6)*pTrunc->local()->rtp()->sTh(th,n);

               //   velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (-6.770275002573076 + 9.0270333367641*r2)*cPh(ph);
               //   velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += (-22.56758334191025 + 90.270333367641*r2 - 72.21626669411279*r4)*sPh(ph);
               //   velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += r*(-7.047567299381681 + 25.691387272599908*r2)*pTrunc->local()->rtp()->sTh(th,n);
               //   velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) *= r*pTrunc->local()->rtp()->cTh(th,n);
               

               /////////////////////////////////////////////////////////////////////////////////////////////////////
               // 
               // FSBenchmark: Thermal convection at E=3e-4, Ra=95, Pr=1, no-slip, fixed temperature
               //

               // Add exact approximation for T10
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (-1.17079 - 2.15284*r2+10.1407*r4-6.89408*r6)*r*sTh;
               // Add coars approximation for T10
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (-1.0- 2.0*r2+10.0*r4-7.0*r6)*r*sTh;

               // Add exact approximation for T43
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = ((-13.056755791851756 +144.99068361538232*r2-131.5648250792794*r4)*c3Ph + (-197.3295714713154 + 360.8499574657975*r2 - 165.31110461974518*r4)*s3Ph)*r4*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((-65.7765238237718+120.28331915526583*r2-55.103701539915065*r4 + (-131.55304764754362+240.5666383105317*r2-110.20740307983012*r4)*c2Th)*c3Ph+(4.3522519306172525-48.33022787179411*r2+43.854941693093124*r4 + (8.704503861234507-96.66045574358822*r2+87.70988338618626*r4)*c2Th)*s3Ph)*r4*sTh*sTh;
               // Add coarse approximation for T43
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = ((-13.0 +145.0*r2-131.5*r4)*c3Ph + (-197.5 + 361.0*r2 - 165.5*r4)*s3Ph)*r4*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((-66.0+120.5*r2-55.0*r4 + (-131.5+240.5*r2-110.0*r4)*c2Th)*c3Ph+(4.5-48.5*r2+44.0*r4 + (8.5-96.5*r2+87.5*r4)*c2Th)*s3Ph)*r4*sTh*sTh;

               // Add exact approximation for P33
               //velV.rOc().rRTP().rR().rShell(n).col(th)(ph) = ((-1.0764226057208708+6.463135000221041*r2-5.537149051028024*r4 + (1.0764226057208703-6.463135000221041*r2+5.537149051028024*r4)*c2Th)*c3Ph + (-20.01647737226051+43.9553626667923*r2-24.070096465811883*r4 + (20.01647737226051-43.9553626667923*r2+24.070096465811886*r4)*c2Th)*s3Ph)*r2*sTh;
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ((-2.152845211441741+19.389405000663125*r2-22.148596204112096*r4)*c3Ph + (-40.032954744521014+131.8660880003769*r2-96.28038586324755*r4)*s3Ph)*r2*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((-40.032954744521014+131.8660880003769*r2-96.28038586324755*r4)*c3Ph + (2.152845211441741-19.389405000663125*r2+22.148596204112096*r4)*s3Ph)*r2*sTh*sTh;
               // Add coarse approximation for P33
               //velV.rOc().rRTP().rR().rShell(n).col(th)(ph) = ((-1.0+6.5*r2-5.5*r4 + (1.0-6.5*r2+5.5*r4)*c2Th)*c3Ph + (-20.0+44.0*r2-24.0*r4 + (20.0-44.0*r2+24.0*r4)*c2Th)*s3Ph)*r2*sTh;
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ((-2.0+19.5*r2-22.0*r4)*c3Ph + (-40.0+132.0*r2-96.5*r4)*s3Ph)*r2*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((-40.0+132.0*r2-96.5*r4)*c3Ph + (2.0-19.5*r2+22.0*r4)*s3Ph)*r2*sTh*sTh;
               //
               /////////////////////////////////////////////////////////////////////////////////////////////////////
               

               /////////////////////////////////////////////////////////////////////////////////////////////////////
               // 
               // FSBenchmark: Thermal convection at E=3e-4, Ra=95, Pr=1, stress-free, fixed temperature
               //

               // Add exact approximation for T10
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (-4.533003573673263 +9.884312649745155*r2-3.281465479238873*r4-1.5395925564120991*r6)*r*sTh;
               // Add coars approximation for T10
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (-4.5 +10.0*r2-3.5*r4-1.5*r6)*r*sTh;
               // Add very coarse approximation for T10
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) = (-5.0 +10.0*r2-3.0*r4-2.0*r6)*r*sTh;
               
               // Add exact approximation for T43
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = ((209.92895118333936 - 461.6697520692621*r2 + 247.7470969880216*r4)*c3Ph + (56.54695912827047- 61.71811873181641*r2 + 8.72871818002742*r4)*s3Ph)*r4*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((18.84898637609016- 20.572706243938804*r2 + 2.9095727266758065*r4 + (37.69797275218031- 41.14541248787761*r2 + 5.819145453351613*r4)*c2Th)*c3Ph + (-69.97631706111312 + 153.8899173564207*r2 - 82.58236566267388*r4 + (-139.95263412222624 + 307.77983471284136*r2 - 165.16473132534773*r4)*c2Th)*s3Ph)*r4*sTh*sTh;
               // Add coarse approximation for T43
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = ((210.0 - 461.5*r2 + 247.5*r4)*c3Ph + (56.5- 61.5*r2 + 8.5*r4)*s3Ph)*r4*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((19.0- 20.5*r2 + 3.0*r4 + (37.5- 41.0*r2 + 6.0*r4)*c2Th)*c3Ph + (-70.0 + 154.0*r2 - 82.5*r4 + (-140.0 + 308.0*r2 - 165.0*r4)*c2Th)*s3Ph)*r4*sTh*sTh;
               // Add very coarse approximation for T43
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) = ((210.0 - 462.0*r2 + 248.0*r4)*c3Ph + (57.0- 62.0*r2 + 9.0*r4)*s3Ph)*r4*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((19.0- 21.0*r2 + 3.0*r4 + (38.0- 41.0*r2 + 6.0*r4)*c2Th)*c3Ph + (-70.0 + 154.0*r2 - 83.0*r4 + (-140.0 + 308.0*r2 - 165.0*r4)*c2Th)*s3Ph)*r4*sTh*sTh;

               // Add exact approximation for P33
               //velV.rOc().rRTP().rR().rShell(n).col(th)(ph) = ((16.724188175612316-41.82620487290966*r2 +25.31016065081357*r4 + (-16.724188175612316 + 41.82620487290966*r2 - 25.31016065081357*r4)*c2Th)*c3Ph + (5.773795540354367- 10.069578281170955*r2 +4.236773196132059*r4 + (-5.773795540354367 + 10.069578281170955*r2 - 4.236773196132059*r4)*c2Th)*s3Ph)*r2*sTh;
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ((33.44837635122463- 125.47861461872898*r2 +101.24064260325429*r4)*c3Ph + (11.547591080708735- 30.208734843512858*r2 + 16.947092784528238*r4)*s3Ph)*r2*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((11.547591080708735-30.208734843512858*r2 + 16.947092784528238*r4)*c3Ph + (-33.44837635122463 + 125.47861461872898*r2 - 101.24064260325429*r4)*s3Ph)*r2*sTh*sTh;
               // Add coarse approximation for P33
               //velV.rOc().rRTP().rR().rShell(n).col(th)(ph) = ((16.5-42.0*r2 +25.5*r4 + (-16.5 + 42.0*r2 - 25.5*r4)*c2Th)*c3Ph + (6.0- 10.0*r2 +4.0*r4 + (-6.0 + 10.0*r2 - 4.0*r4)*c2Th)*s3Ph)*r2*sTh;
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ((33.5- 125.5*r2 +101.0*r4)*c3Ph + (11.5- 30.0*r2 + 17.0*r4)*s3Ph)*r2*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((11.5-30.0*r2 + 17.0*r4)*c3Ph + (-33.5 + 125.5*r2 - 101.0*r4)*s3Ph)*r2*sTh*sTh;
               // Add very coarse approximation for P33
               //velV.rOc().rRTP().rR().rShell(n).col(th)(ph) = ((17.0-42.0*r2 +25.0*r4 + (-17.0 + 42.0*r2 - 25.0*r4)*c2Th)*c3Ph + (6.0- 10.0*r2 +4.0*r4 + (-6.0 + 10.0*r2 - 4.0*r4)*c2Th)*s3Ph)*r2*sTh;
               //velV.rOc().rRTP().rTheta().rShell(n).col(th)(ph) += ((33.0- 125.0*r2 +101.0*r4)*c3Ph + (12.0- 30.0*r2 + 17.0*r4)*s3Ph)*r2*cTh*sTh*sTh;
               //velV.rOc().rRTP().rPhi().rShell(n).col(th)(ph) += ((12.0-30.0*r2 + 17.0*r4)*c3Ph + (-33.0 + 125.0*r2 - 101.0*r4)*s3Ph)*r2*sTh*sTh;
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
      velV.rOc().rPerturbation().rTor().rLShell(2)(2,1) = ampT*EPMComplex(-9.94764,0.77007);
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
