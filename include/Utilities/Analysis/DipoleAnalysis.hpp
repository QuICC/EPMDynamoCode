/** \file DipoleAnalysis.hpp
 *  \brief Defines possible analysis to be done on state field
 */

#ifndef DIPOLEANALYSIS_HPP
#define DIPOLEANALYSIS_HPP

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
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/DipoleFile.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits for libration analysis
    */
   class DipoleTrits
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
         static const bool UseSpecMagnetic = true;

         /// Requires Magnetic curl computations
         static const bool UseMagneticCurl = false;

         /// Requires RTP Velocity computations
         static const bool UseRTPVelocity = false;

         /// Requires Velocity computations
         static const bool UseSpecVelocity = false;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;
   };

   /**
    * \brief Defines possible analysis to be done on state field
    */
   template <typename TSimTraits> class DipoleAnalysis: public GeneratorBase<TSimTraits>
   {
      public:
         /// Typdef for the AnalysisTraits type
         typedef DipoleTrits  AnalysisTraits;

         /**
          * @brief Private constructor
          */
         DipoleAnalysis();

         /**
          * @brief Destructor
          */
         virtual ~DipoleAnalysis() {};

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

   template <typename TSimTraits> DipoleAnalysis<TSimTraits>::DipoleAnalysis()
      : GeneratorBase<TSimTraits>()
   {
   }

   template <typename TSimTraits> void DipoleAnalysis<TSimTraits>::analyseRTPCodensity()
   {
      SmartTruncation pTrunc = this->codC().oc().trunc();

      std::cout << "------ Codensity RTP field analysis ------" << std::endl;
   }

   template <typename TSimTraits> void DipoleAnalysis<TSimTraits>::analyseRTPMagnetic()
   {
      SmartTruncation pTrunc = this->magB().oc().trunc();

      std::cout << "------ Magnetic RTP field analysis ------" << std::endl;
   }

   template <typename TSimTraits> void DipoleAnalysis<TSimTraits>::analyseRTPVelocity()
   {
      SmartTruncation pTrunc = this->velV().oc().trunc();

      std::cout << "------ Velocity RTP field analysis ------" << std::endl;
   }

   template <typename TSimTraits> void DipoleAnalysis<TSimTraits>::analyseSpecCodensity()
   {
      SmartTruncation pTrunc = this->codC().oc().trunc();

      std::cout << "------ Codensity spectral field analysis ------" << std::endl;
   }

   template <typename TSimTraits> void DipoleAnalysis<TSimTraits>::analyseSpecMagnetic()
   {
      SmartTruncation pTrunc = this->magB().oc().trunc();

      std::cout << "------ Magnetic spectral field analysis ------" << std::endl;

      SmartASCIIWriter pOutFile;
      
      pOutFile.reset(new DipoleFile<typename TSimTraits::MagType>(this->magB(), "mag", this->mTSParams, this->mTransform.radBasis(), 0));

      pOutFile->init();
      pOutFile->write();
      pOutFile->finalise();

      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<typename TSimTraits::MagType> > pMagSpectrum(new SpectrumFile<typename TSimTraits::MagType>(this->magB(), "mag", 1));

      pMagSpectrum->init();
      pMagSpectrum->write();
      pMagSpectrum->finalise();
   }

   template <typename TSimTraits> void DipoleAnalysis<TSimTraits>::analyseSpecVelocity()
   {
      std::cout << "------ Velocity spectral field analysis ------" << std::endl;
   }
}

#endif // DIPOLEANALYSIS_HPP
