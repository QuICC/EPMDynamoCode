/** \file BasicImposed.hpp
 *  \brief Sources definition function for Sakuraba's paper
 */

#ifndef BASICIMPOSED_HPP
#define BASICIMPOSED_HPP

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
    * @brief Traits of basic imposed generator
    */
   class BasicImposedTraits
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
    * \brief Imposed field definition function
    */
   template <typename TSimTraits> class BasicImposed
   {
      public:
         /// Typdef for the basic imposed traits type
         typedef BasicImposedTraits  ImposedTraits;

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
         BasicImposed();

         /**
          * @brief Destructor
          */
         virtual ~BasicImposed() {};
   };

   template <typename TSimTraits> void BasicImposed<TSimTraits>::setRTPCodensity(typename BasicImposed<TSimTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TSimTraits> void BasicImposed<TSimTraits>::setRTPMagnetic(typename BasicImposed<TSimTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
      }
   }

   template <typename TSimTraits> void BasicImposed<TSimTraits>::setRTPVelocity(typename BasicImposed<TSimTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
      }
   }

   template <typename TSimTraits> void BasicImposed<TSimTraits>::setSpecCodensity(typename BasicImposed<TSimTraits>::Codensity &codC)
   {
   }

   template <typename TSimTraits> void BasicImposed<TSimTraits>::setSpecMagnetic(typename BasicImposed<TSimTraits>::Magnetic &magB)
   {
   }

   template <typename TSimTraits> void BasicImposed<TSimTraits>::setSpecVelocity(typename BasicImposed<TSimTraits>::Velocity &velV)
   {
   }

   template <typename TSimTraits> BasicImposed<TSimTraits>::BasicImposed()
   {
   }

}

#endif // BASICIMPOSED_HPP
