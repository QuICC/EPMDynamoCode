/** \file SerialTimer.hpp
 *  \brief Implementation of a serial timer
 */

#ifndef SERIALTIMER_HPP
#define SERIALTIMER_HPP

// System includes
//
#include <time.h>
#if defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

// External includes
//

// Project includes
//
#include "General/TimerBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a serial timer
    */
   class SerialTimer: public TimerBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         SerialTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         virtual ~SerialTimer() {};

         /**
          * @brief Start clock
          */
         virtual void start();

         /**
          * @brief Stop clock
          */
         virtual void stop();

         /**
          * @brief Get elapsed time
          */
         virtual EPMFloat time() const;

         /**
          * @brief Reset timer (stop and restart)
          */
         virtual EPMFloat resetTimer();
         
      protected:

      private:
         /**
          * @brief Start
          */
         timespec mStart;

         /**
          * @brief Stop
          */
         timespec mStop;

         /**
          * @brief Compute elapsed seconds between start and stop
          */
         EPMFloat elapsedSeconds() const;
   };

   /// Typedef for a generic Timer
   typedef SerialTimer  TimerType;

}

#endif // SERIALTIMER_HPP
