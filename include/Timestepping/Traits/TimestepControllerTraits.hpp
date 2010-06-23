/** \file TimestepControllerTraits.hpp
 *  \brief Timestep controller traits for the predefined controllers
 */

#ifndef TIMESTEPCONTROLLERTRAITS_HPP
#define TIMESTEPCONTROLLERTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Timestep controller traits
    *
    * \tparam TType The type of timestep controller
    */
   template <TimestepCtrlTypes TType> class TimestepControllerTraits;

   template <> class TimestepControllerTraits<ElementaryCtrl>
   {
      public:
         /**
          * @brief \f$\beta_1\f$ parameter of the elementary controller
          */
         static const EPMFloat KBETA1;

         /**
          * @brief \f$\beta_2\f$ parameter of the elementary controller
          */
         static const EPMFloat KBETA2;

         /**
          * @brief \f$\alpha\f$ parameter of the elementary controller
          */
         static const EPMFloat ALPHA;
   };

   template <> class TimestepControllerTraits<PI42Ctrl>
   {
      public:
         /**
          * @brief \f$\beta_1\f$ parameter of the PI4.2 controller
          */
         static const EPMFloat KBETA1;

         /**
          * @brief \f$\beta_2\f$ pf the PI4.2 controller
          */
         static const EPMFloat KBETA2;

         /**
          * @brief \f$\alpha\f$ pa the PI4.2 controller
          */
         static const EPMFloat ALPHA;
   };

   template <> class TimestepControllerTraits<H211BCtrl>
   {
      public:
         /**
          * @brief \f$\beta_1\f$ pf the H211b controller
          */
         static const EPMFloat KBETA1;

         /**
          * @brief \f$\beta_2\f$ pf the H211b controller
          */
         static const EPMFloat KBETA2;

         /**
          * @brief \f$\alpha\f$ pa the H211b controller
          */
         static const EPMFloat ALPHA;
   };

}

#endif // TIMESTEPCONTROLLERTRAITS_HPP
