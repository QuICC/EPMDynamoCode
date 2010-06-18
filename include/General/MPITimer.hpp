/** \file MPITimer.hpp
 *  \brief Implementation of a MPI timer
 */

#ifndef MPITIMER_HPP
#define MPITIMER_HPP

// System includes
//
#include <mpi.h>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/TimerBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a MPI timer
    */
   class MPITimer: public TimerBase
   {
      public:
         /**
          * @brief Constructor
          */
         MPITimer();

         /**
          * @brief Destructor
          */
         virtual ~MPITimer() {};

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
         EPMFloat mStart;

         /**
          * @brief Stop
          */
         EPMFloat mStop;
   };

   /// Typedef for a generic Timer
   typedef MPITimer  TimerType;

}

#endif // MPITIMER_HPP
