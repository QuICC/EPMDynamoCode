/** \file LibrationAnalysis.hpp
 *  \brief Defines possible analysis to be done on state field
 */

#ifndef LIBRATIONANALYSIS_HPP
#define LIBRATIONANALYSIS_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "Utilities/GeneratorBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits for libration analysis
    */
   class LibrationTraits
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
         static const bool UseRTPVelocity = false;

         /// Requires Velocity computations
         static const bool UseSpecVelocity = true;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;
   };

   /**
    * \brief Defines possible analysis to be done on state field
    */
   template <typename TSimTraits> class LibrationAnalysis: public GeneratorBase<TSimTraits>
   {
      public:
         /// Typdef for the AnalysisTraits type
         typedef LibrationTraits  AnalysisTraits;

         /**
          * @brief Private constructor
          */
         LibrationAnalysis();

         /**
          * @brief Destructor
          */
         virtual ~LibrationAnalysis() {};

         /**
          * @brief Analyse RTP value for the codensity scalar
          */
         void analyseRTPCodensity();

         /**
          * @brief Analyse RTP value for the magnetic field
          */
         void analyseRTPMagnetic();

         /**
          * @brief Analyse RTP value for the velocity field
          */
         void analyseRTPVelocity();

         /**
          * @brief Analyse spectral coefficients for the codensity scalar
          */
         void analyseSpecCodensity();

         /**
          * @brief Analyse spectral coefficients for the velocity field
          */
         void analyseSpecMagnetic();

         /**
          * @brief Analyse spectral coefficients for the magnetic field
          */
         void analyseSpecVelocity();

      private:
   };

   template <typename TSimTraits> LibrationAnalysis<TSimTraits>::LibrationAnalysis()
      : GeneratorBase<TSimTraits>()
   {
   }

   template <typename TSimTraits> void LibrationAnalysis<TSimTraits>::analyseRTPCodensity()
   {
      SmartTruncation pTrunc = this->codC().oc().trunc();

      std::cout << "------ Codensity RTP field analysis ------" << std::endl;
   }

   template <typename TSimTraits> void LibrationAnalysis<TSimTraits>::analyseRTPMagnetic()
   {
      SmartTruncation pTrunc = this->magB().oc().trunc();

      std::cout << "------ Magnetic RTP field analysis ------" << std::endl;
   }

   template <typename TSimTraits> void LibrationAnalysis<TSimTraits>::analyseRTPVelocity()
   {
      SmartTruncation pTrunc = this->velV().oc().trunc();

      std::cout << "------ Velocity RTP field analysis ------" << std::endl;
   }

   template <typename TSimTraits> void LibrationAnalysis<TSimTraits>::analyseSpecCodensity()
   {
      SmartTruncation pTrunc = this->codC().oc().trunc();

      std::cout << "------ Codensity spectral field analysis ------" << std::endl;

      // Update and output the energy spectra
      this->codC().rOc().updateSpectra();
      std::cout << "\t" << "Energy: " << this->codC().oc().energy().transpose() << std::endl;
   }

   template <typename TSimTraits> void LibrationAnalysis<TSimTraits>::analyseSpecMagnetic()
   {
      SmartTruncation pTrunc = this->magB().oc().trunc();

      std::cout << "------ Magnetic spectral field analysis ------" << std::endl;

      // Update the energy spectra
      this->magB().rOc().updateSpectra();
      std::cout <<  "\t" << "Energy: " << this->magB().oc().energy().transpose() << std::endl;
   }

   template <typename TSimTraits> void LibrationAnalysis<TSimTraits>::analyseSpecVelocity()
   {
      SmartTruncation pTrunc = this->velV().oc().trunc();

      std::cout << "------ Velocity spectral field analysis ------" << std::endl;

      // Update the energy spectra
      this->velV().rOc().updateSpectra();
      std::cout <<  "\t" << "Energy: " << this->velV().oc().energy().transpose() << std::endl;

      // Compute scalar product with solid body rotations
      std::cout <<  "\t" << "Projection on solid body rotations " << std::endl;

      // Compute scalar product with solid body rotation around x axis
      ArrayZ xAxis = ArrayZ::Zero(pTrunc->sim()->rad()->nN());
      xAxis(0).real() = -0.443113;
      EPMFloat xProj = 0.0;

      // Compute scalar product with solid body rotation around y axis
      ArrayZ yAxis = ArrayZ::Zero(pTrunc->sim()->rad()->nN());
      yAxis(0).imag() = 0.443113;
      EPMFloat yProj = 0.0;

      int l = 1;
      int l_ = 1;
      int m = 1;
      EPMFloat shFactor = 4.0*4.0 * MathConstants::PI / static_cast<EPMFloat>(2*l_+1);
      EPMFloat lfactor = static_cast<EPMFloat>(l_*(l_+1));
      int nN = pTrunc->sim()->rad()->nN();

      for(int n = 0; n < nN; ++n)
      {
         for(int k = 0; k < nN; ++k)
         {
            xProj += this->mTransform.radBasis().at(l).eWeights()(k,n) * (this->velV().oc().perturbation().tor().lshell(l)(n,m).real()*xAxis(k).real() + this->velV().oc().perturbation().tor().lshell(l)(n,m).imag()*xAxis(k).imag());
            yProj += this->mTransform.radBasis().at(l).eWeights()(k,n) * (this->velV().oc().perturbation().tor().lshell(l)(n,m).real()*yAxis(k).real() + this->velV().oc().perturbation().tor().lshell(l)(n,m).imag()*yAxis(k).imag());
         }
      }

      xProj *= shFactor*lfactor;
      yProj *= shFactor*lfactor;

      std::cout <<  "\t" << "x axis: " << xProj << std::endl;
      std::cout <<  "\t" << "y axis: " << yProj << std::endl;
   }
}

#endif // LIBRATIONANALYSIS_HPP
