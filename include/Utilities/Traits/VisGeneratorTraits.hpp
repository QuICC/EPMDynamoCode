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
         static const bool UseSpecCodensity = TVisTraits::VisCodensity || TVisTraits::VisCodensityGrad;

         static const bool UseCodensityGrad = TVisTraits::VisCodensityGrad;


         static const bool UseMagneticCurl = TVisTraits::VisMagnetic || TVisTraits::VisLorentz;

         static const bool UseSpecMagnetic = TVisTraits::VisMagnetic;


         static const bool UseSpecVelocity = TVisTraits::VisVelocity || TVisTraits::VisVorticity || TVisTraits::VisHelicity;

         static const bool UseVelocityCurl = TVisTraits::VisVorticity || TVisTraits::VisHelicity;



         static const bool UseRTPCodensity = false;

         static const bool UseRTPMagnetic = false;

         static const bool UseRTPVelocity = false;
   };

}

#endif // VISGENERATORTRAITS_HPP
