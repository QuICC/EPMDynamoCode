/** \file SerialTimer.cpp
 *  \brief Source of the serial timer implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "General/SerialTimer.hpp"

// Project includes
//

namespace EPMDynamo {

   SerialTimer::SerialTimer(const bool autostart)
   {
      // Check if timer should be started
      if(autostart)
      {
         this->start();
      }
   }

   void SerialTimer::start()
   {
      #ifdef __APPLE__
      mach_absolute_time();
      #else
      clock_gettime(CLOCK_REALTIME, &this->mStart);
      #endif
   }

   void SerialTimer::stop()
   {
       #ifdef __APPLE__
       mach_absolute_time();
       #else
       clock_gettime(CLOCK_REALTIME, &this->mStop);
       #endif
   }

   EPMFloat SerialTimer::time() const
   {
      return elapsedSeconds();
   }

   EPMFloat SerialTimer::resetTimer()
   {
      this->stop();
      EPMFloat tmp = this->time();
      this->mStart = this->mStop;

      return tmp;
   }

   EPMFloat SerialTimer::elapsedSeconds() const
   {
      return static_cast<EPMFloat>(this->mStop.tv_sec - this->mStart.tv_sec) + static_cast<EPMFloat>(this->mStop.tv_nsec - this->mStart.tv_nsec)/1.0e9;
   }

}
