/** \file InfluenceTTraits.hpp
 *  \brief Timestepper traits including influence matrix step(s)
 */

#ifndef INFLUENCETTRAITS_HPP
#define INFLUENCETTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Timestepper traits including influence matrix step(s)
    *
    * \tparam TTimeTraits Timestepping traits
    */
   template <typename TTimeTraits> class InfluenceTTraits
   {
      public:
         /// Define typedef to use the influence matrix step in timestepper
         typedef typename TTimeTraits::InfluenceTimestepper   Timestepper;  
   };

}

#endif // INFLUENCETTRAITS_HPP
