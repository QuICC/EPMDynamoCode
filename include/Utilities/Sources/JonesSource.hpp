/** \file JonesSource.hpp
 *  \brief Sources definition function for Jones's paper
 */

#ifndef JONESSOURCE_HPP
#define JONESSOURCE_HPP

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
    * \brief Sources definition function for Jones's paper
    */
   template <typename TGenTraits> class JonesSource
   {
      public:
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
         JonesSource();

         /**
          * @brief Destructor
          */
         virtual ~JonesSource() {};
   };

   template <typename TGenTraits> void JonesSource<TGenTraits>::setRTPCodensity(typename JonesSource<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void JonesSource<TGenTraits>::setRTPMagnetic(typename JonesSource<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void JonesSource<TGenTraits>::setRTPVelocity(typename JonesSource<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
      }
   }

   template <typename TGenTraits> void JonesSource<TGenTraits>::setSpecCodensity(typename JonesSource<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      codC.rOcSrc().rLShell(0)(0,0) = 3.0;
   }

   template <typename TGenTraits> void JonesSource<TGenTraits>::setSpecMagnetic(typename JonesSource<TGenTraits>::Magnetic &magB)
   {
   }

   template <typename TGenTraits> void JonesSource<TGenTraits>::setSpecVelocity(typename JonesSource<TGenTraits>::Velocity &velV)
   {
   }

   template <typename TGenTraits> JonesSource<TGenTraits>::JonesSource()
   {
   }

}

#endif // JONESSOURCE_HPP
