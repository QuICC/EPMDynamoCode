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

   SerialTimer::SerialTimer()
   {
   }

   void SerialTimer::start()
   {
      clock_gettime(CLOCK_REALTIME, &mStart);
   }

   void SerialTimer::stop()
   {
      clock_gettime(CLOCK_REALTIME, &mStop);
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
