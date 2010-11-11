/** \file VisGeneratorTraits.cpp
 *  \brief Convert Visualisation traits to generator traits
 */

#ifndef VISGENERATORTRAITS_HPP
#define VISGENERATORTRAITS_HPP

namespace EPMDynamo {

   /**
    * \brief Convert Visualisation traits to generator traits
    *
    * \tparam TVisTraits Traits of visualisation output
    */
   template <typename TVisTraits>  class VisGeneratorTraits
   {
      public:
         /// Computation of the RTP transform for the codensity scalar is required?
         static const bool UseSpecCodensity = TVisTraits::VisCodensity || TVisTraits::VisCodensityGrad;

         /// Computation of the RTP transform for the codensity gradient is required?
         static const bool UseCodensityGrad = TVisTraits::VisCodensityGrad;


         /// Computation of the RTP transform for the magnetic curl is required?
         static const bool UseMagneticCurl = TVisTraits::VisMagnetic || TVisTraits::VisLorentz;

         /// Computation of the RTP transform for the magnetic field is required?
         static const bool UseSpecMagnetic = TVisTraits::VisMagnetic;


         /// Computation of the RTP transform for the velocity field is required?
         static const bool UseSpecVelocity = TVisTraits::VisVelocity || TVisTraits::VisVorticity || TVisTraits::VisHelicity;

         /// Computation of the RTP transform for the velocity curl is required?
         static const bool UseVelocityCurl = TVisTraits::VisVorticity || TVisTraits::VisHelicity;



         /// Computation of the spectral transform for the codensity scalar is required? (is mostly never going to be true)
         static const bool UseRTPCodensity = false;

         /// Computation of the spectral transform for the magnetic field is required? (is mostly never going to be true)
         static const bool UseRTPMagnetic = false;

         /// Computation of the spectral transform for the velocity field is required? (is mostly never going to be true)
         static const bool UseRTPVelocity = false;
   };

}

#endif // VISGENERATORTRAITS_HPP
