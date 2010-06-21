/** \file SimpleTTraits.hpp
 *  \brief Simple Timestepper traits without influence matrix step(s)
 */

#ifndef SIMPLETTRAITS_HPP
#define SIMPLETTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Simple Timestepper traits without influence matrix step(s)
    *
    * \tparam TTimeTraits Timestepping traits
    */
   template <typename TTimeTraits> class SimpleTTraits
   {
      public:
         /// Define typedef for standard timestepper
         typedef typename TTimeTraits::Timestepper   Timestepper;  
   };

}

#endif // SIMPLETTRAITS_HPP
