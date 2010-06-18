/** \file MPITimer.cpp
 *  \brief Source of the implementation of MPI timer
 */

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// External includes
//

// Class include
//
#include "General/MPITimer.hpp"

// Project includes
//

namespace EPMDynamo {

   MPITimer::MPITimer()
      : mStart(0.0), mStop(0.0)
   {
   }

   void MPITimer::start()
   {
      this->mStart = MPI_Wtime();
   }

   void MPITimer::stop()
   {
      this->mStop = MPI_Wtime();
   }

   EPMFloat MPITimer::time() const
   {
      return this->mStop - this->mStart;
   }

   EPMFloat MPITimer::resetTimer()
   {
      this->stop();
      EPMFloat tmp = this->time();
      this->mStart = this->mStop;

      return tmp;
   }

}
