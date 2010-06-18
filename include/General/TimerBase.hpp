/** \file TimerBase.hpp
 *  \brief Base of the implementation of a timer
 */

#ifndef TIMERBASE_HPP
#define TIMERBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a timer
    */
   class TimerBase
   {
      public:
         /**
          * @brief Constructor
          */
         TimerBase();

         /**
          * @brief Destructor
          */
         virtual ~TimerBase() {};

         /**
          * @brief Start clock
          */
         virtual void start() = 0;

         /**
          * @brief Stop clock
          */
         virtual void stop() = 0;

         /**
          * @brief Get elapsed time
          */
         virtual EPMFloat time() const = 0;

         /**
          * @brief Reset timer (stop and restart)
          */
         virtual EPMFloat resetTimer() = 0;
         
      protected:

      private:
   };

}

#endif // TIMERBASE_HPP
