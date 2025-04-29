/** \file ExtractRadiusAnalysis.hpp
 *  \brief Defines possible analysis to be done on state field
 */

#ifndef EXTRACTRADIUSANALYSIS_HPP
#define EXTRACTRADIUSANALYSIS_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iomanip>

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/MathConstants.hpp"
#include "Domain/Truncation.hpp"
#include "Utilities/GeneratorBase.hpp"
#include "IO/ASCII/ExtractRadiusFile.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits for libration analysis
    */
   class ExtractRadiusTraits
   {
      public:
         /// Requires RTP Codensity computations
         static const bool UseRTPCodensity = true;

         /// Requires spectral Codensity computations
         static const bool UseSpecCodensity = true;

         /// Requires Codensity gradient computations
         static const bool UseCodensityGrad = false;

         /// Requires RTP Magnetic computations
         static const bool UseRTPMagnetic = true;

         /// Requires spectral Magnetic computations
         static const bool UseSpecMagnetic = true;

         /// Requires Magnetic curl computations
         static const bool UseMagneticCurl = false;

         /// Requires RTP Velocity computations
         static const bool UseRTPVelocity = true;

         /// Requires Velocity computations
         static const bool UseSpecVelocity = true;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;

         /// Colatitude at which to extract the radius (degrees)
         static const int colatitude = 30;

         /// Longitude at which to extract the radius (degrees)
         static const int longitude = -1;
   };

   /**
    * \brief Defines possible analysis to be done on state field
    */
   template <typename TSimTraits> class ExtractRadiusAnalysis: public GeneratorBase<TSimTraits>
   {
      public:
         /// Typdef for the AnalysisTraits type
         typedef ExtractRadiusTraits  AnalysisTraits;

         /**
          * @brief Private constructor
          */
         ExtractRadiusAnalysis();

         /**
          * @brief Destructor
          */
         virtual ~ExtractRadiusAnalysis() {};

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
         /**
          * @brief Theta point to extract
          */
         int mTheta;

         /**
          * @brief Phi point to extract
          */
         int mPhi;

         /**
          * @brief Setup
          */
         void setup(SmartTruncation pTrunc);
   };

   template <typename TSimTraits> ExtractRadiusAnalysis<TSimTraits>::ExtractRadiusAnalysis()
      : GeneratorBase<TSimTraits>()
   {

   }

   template <typename TSimTraits> void ExtractRadiusAnalysis<TSimTraits>::setup(SmartTruncation pTrunc)
   {
      std::cout << "----------------- Radius extraction ------------------" << std::endl;
      std::cout << "   -> Preparing extraction of values along a radius   " << std::endl;
      for(int i = 0; i < pTrunc->sim()->hoz()->nTh(); i++)
      {
         if(ExtractRadiusTraits::colatitude/180.*MathConstants::PI > pTrunc->sim()->hoz()->thGrid()(i))
         {
            this->mTheta = i;
            break;
         }
      }

      if(ExtractRadiusTraits::longitude < 0)
      {
         this->mPhi = -1;
      } else
      {
         for(int i = 0; i < pTrunc->sim()->hoz()->nPh(); i++)
         {
            if(ExtractRadiusTraits::longitude/360.*MathConstants::PI < pTrunc->sim()->hoz()->phGrid()(i))
            {
               this->mPhi = i;
               break;
            }
         }
      }

      std::cout << "   -> Selected latitude: " << (pTrunc->sim()->hoz()->thGrid()(this->mTheta)/MathConstants::PI)*180. << "(" << this->mTheta << ")" << std::endl;
      if(this->mPhi < 0)
      {
         std::cout << "   -> All longitudes are being extracted" << std::endl;
      } else
      {
         std::cout << "   -> Selected longitude: " << (pTrunc->sim()->hoz()->phGrid()(this->mPhi)/MathConstants::PI)*180. << "(" << this->mPhi << ")" << std::endl;
      }
   }

   template <typename TSimTraits> void ExtractRadiusAnalysis<TSimTraits>::analyseRTPCodensity()
   {
      SmartTruncation pTrunc = this->codC().oc().trunc();

      this->setup(pTrunc);

//      SmartASCIIWriter pOutFile(new ExtractRadiusFile<typename TSimTraits::CodType>(this->codC(), this->mTheta, this->mPhi, "cod", 0));

//      pOutFile->init();
//      pOutFile->write();
//      pOutFile->finalise();
   }

   template <typename TSimTraits> void ExtractRadiusAnalysis<TSimTraits>::analyseRTPMagnetic()
   {
      SmartTruncation pTrunc = this->magB().oc().trunc();

      this->setup(pTrunc);

      SmartASCIIWriter pOutFile;
      
      pOutFile.reset(new ExtractRadiusFile<typename TSimTraits::MagType>(this->magB(), this->mTheta, this->mPhi, "magR", 0));

      pOutFile->init();
      pOutFile->write();
      pOutFile->finalise();
      
      pOutFile.reset(new ExtractRadiusFile<typename TSimTraits::MagType>(this->magB(), this->mTheta, this->mPhi, "magTh", 1));

      pOutFile->init();
      pOutFile->write();
      pOutFile->finalise();
      
      pOutFile.reset(new ExtractRadiusFile<typename TSimTraits::MagType>(this->magB(), this->mTheta, this->mPhi, "magPh", 2));

      pOutFile->init();
      pOutFile->write();
      pOutFile->finalise();
   }

   template <typename TSimTraits> void ExtractRadiusAnalysis<TSimTraits>::analyseRTPVelocity()
   {
      SmartTruncation pTrunc = this->velV().oc().trunc();

      this->setup(pTrunc);

      SmartASCIIWriter pOutFile;
      
      pOutFile.reset(new ExtractRadiusFile<typename TSimTraits::VelType>(this->velV(), this->mTheta, this->mPhi, "velR", 0));

      pOutFile->init();
      pOutFile->write();
      pOutFile->finalise();
      
      pOutFile.reset(new ExtractRadiusFile<typename TSimTraits::VelType>(this->velV(), this->mTheta, this->mPhi, "velTh", 1));

      pOutFile->init();
      pOutFile->write();
      pOutFile->finalise();
      
      pOutFile.reset(new ExtractRadiusFile<typename TSimTraits::VelType>(this->velV(), this->mTheta, this->mPhi, "velPh", 2));

      pOutFile->init();
      pOutFile->write();
      pOutFile->finalise();
   }

   template <typename TSimTraits> void ExtractRadiusAnalysis<TSimTraits>::analyseSpecCodensity()
   {
   }

   template <typename TSimTraits> void ExtractRadiusAnalysis<TSimTraits>::analyseSpecMagnetic()
   {
   }

   template <typename TSimTraits> void ExtractRadiusAnalysis<TSimTraits>::analyseSpecVelocity()
   {
   }
}

#endif // EXTRACTRADIUSANALYSIS_HPP
