/** \file ExecutionTimer.hpp
 *  \brief Implementation of a execution timer
 */

#ifndef EXECUTIONTIMER_HPP
#define EXECUTIONTIMER_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

#ifdef EPMDYNAMO_MPI
   #include "General/MPITimer.hpp"
#else
   #include "General/SerialTimer.hpp"
#endif // EPMDYNAMO_MPI

namespace EPMDynamo {

   /**
    * \brief Implementation of a execution timer
    */
   class ExecutionTimer: public TimerType
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         ExecutionTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         virtual ~ExecutionTimer() {};

         /**
          * @brief Update the timing of given region
          *
          * @param id Id of the timed region
          */
         void update(const int id);

         /**
          * @brief Analyze the measured times
          */
         void analyze();

         /**
          * @brief Print execution time information
          */
         void printInfo() const;
         
      protected:

      private:
         /**
          * @brief Mean execution times
          */
         Array mMeanTimes;

         /**
          * @brief Max execution times
          */
         Array mMaxTimes;

         /**
          * @brief Mean execution times
          */
         Array mMinTimes;
   };

}

#endif // EXECUTIONTIMER_HPP
