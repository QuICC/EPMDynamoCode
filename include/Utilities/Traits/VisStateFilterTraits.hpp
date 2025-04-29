/** \file VisStateFilterTraits.cpp
 *  \brief Convert Visualisation traits to state filter traits
 */

#ifndef VISSTATEFILTERTRAITS_HPP
#define VISSTATEFILTERTRAITS_HPP

namespace EPMDynamo {

   /**
    * \brief Convert Visualisation traits to state filter traits
    *
    * \tparam TVisTraits Traits of visualisation output
    */
   template <typename TVisTraits, const int Filter>  class VisStateFilterTraits
   {
      public:
         /// Specify part of codensity to read
         static const int ReadCodensity = TVisTraits::VisCodensity;

         /// Specify part of magnetic to read
         static const int ReadMagnetic = Filter*TVisTraits::VisMagnetic;

         /// Specify part of velocity field to read
         static const int ReadVelocity = Filter*TVisTraits::VisVelocity;
   };

}

#endif // VISSTATEFILTERTRAITS_HPP
