/** \file SakurabaSource.hpp
 *  \brief Sources definition function for Sakuraba's paper
 */

#ifndef SAKURABASOURCE_HPP
#define SAKURABASOURCE_HPP

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
    * @brief Traits of Sakuraba source generator
    */
   class SakurabaSourceTraits
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
         static const bool UseSpecVelocity = false;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;
   };

   /**
    * \brief Sources definition function for Sakuraba's paper
    */
   template <typename TSimTraits> class SakurabaSource
   {
      public:
         /// Typdef for the Sakuraba's paper source traits type
         typedef SakurabaSourceTraits  SourceTraits;

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
         SakurabaSource();

         /**
          * @brief Destructor
          */
         virtual ~SakurabaSource() {};
   };

   template <typename TSimTraits> void SakurabaSource<TSimTraits>::setRTPCodensity(typename SakurabaSource<TSimTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TSimTraits> void SakurabaSource<TSimTraits>::setRTPMagnetic(typename SakurabaSource<TSimTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
      }
   }

   template <typename TSimTraits> void SakurabaSource<TSimTraits>::setRTPVelocity(typename SakurabaSource<TSimTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
      }
   }

   template <typename TSimTraits> void SakurabaSource<TSimTraits>::setSpecCodensity(typename SakurabaSource<TSimTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      codC.rOcSrc().rLShell(0)(0,0) = 3.0*10.0;
   }

   template <typename TSimTraits> void SakurabaSource<TSimTraits>::setSpecMagnetic(typename SakurabaSource<TSimTraits>::Magnetic &magB)
   {
   }

   template <typename TSimTraits> void SakurabaSource<TSimTraits>::setSpecVelocity(typename SakurabaSource<TSimTraits>::Velocity &velV)
   {
   }

   template <typename TSimTraits> SakurabaSource<TSimTraits>::SakurabaSource()
   {
   }

}

#endif // SAKURABASOURCE_HPP
