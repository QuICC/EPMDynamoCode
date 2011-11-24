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
      std::cout << "------ Velocity spectral field analysis ------" << std::endl;

      // Update the energy spectra
      this->velV().rOc().updateSpectra();
      std::cout <<  "\t" << "Energy: " << this->velV().oc().energy().transpose() << std::endl;

      // Compute scalar product with solid body rotations
      std::cout <<  "\t" << "Projection on solid body rotations " << std::endl;

      Array solidProjection;

      solidProjection = this->velV().rOc().rPerturbation().computeXYZSolidProjection(this->mTransform.radBasis());

      std::cout <<  "\t" << "x axis: " << solidProjection(0) << std::endl;
      std::cout <<  "\t" << "y axis: " << solidProjection(1) << std::endl;
      std::cout <<  "\t" << "z axis: " << solidProjection(2) << std::endl;
      std::cout <<  "\t" << "norm: " << solidProjection(3) << std::endl;
      std::cout <<  "\t" << "poincare: " << std::sqrt(solidProjection(0)*solidProjection(0) + solidProjection(1)*solidProjection(1)) << std::endl;
   }
}

#endif // LIBRATIONANALYSIS_HPP
