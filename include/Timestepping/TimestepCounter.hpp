/** \file TimestepCounter.hpp
 *  \brief Implementation of a timestep related data counter
 */

#ifndef TIMESTEPCOUNTER_HPP
#define TIMESTEPCOUNTER_HPP

// System includes
//
#include <iostream>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a timestep related data counter
    */
   class TimestepCounter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param maxtstep Maximum number of timesteps
          * @param asciiRate ASCII file save rate
          * @param stateRate State file save rate
          */
         TimestepCounter(int maxtstep, int asciiRate, int stateRate);

         /**
          * @brief Destructor
          */
         virtual ~TimestepCounter() {};

         /**
          * @brief keep running simulation?
          */
         bool keepRunning() const;

         /**
          * @brief Increment counter
          *
          * @param time Time for time trigger
          */
         void increment(const EPMFloat time);

         /**
          * @brief Get the number of steps
          */
         int steps() const;

         /**
          * @brief is maxtstep unlimited?
          */
         bool isUnlimited() const;

         /**
          * @brief Trigger ASCII save?
          */
         bool triggerASCII() const;

         /**
          * @brief Trigger state save?
          */
         bool triggerState() const;

         /**
          * @brief Print some information about run
          */
         void printInfo() const;
         
      protected:

      private:

         /**
          * @brief Is IO step triggered or time for ASCII?
          */
         const bool mUseASCIIStep;

         /**
          * @brief Is IO step triggered or time for State?
          */
         const bool mUseStateStep;

         /**
          * @brief Maxtstep
          */
         const int   mMaxtstep;

         /**
          * @brief Step counter
          */
         int   mSteps;

         /**
          * @brief Save rate for ASCII files
          */
         int   mARate;

         /**
          * @brief Save rate for State files
          */
         int   mSRate;

         /**
          * @brief Has ASCII IO been triggered?
          */
         bool mIsATriggered;

         /**
          * @brief Has ASCII IO been triggered?
          */
         bool mIsSTriggered;

         /**
          * @brief Intermediate data for time trigger for ASCII
          */
         int mNextAWrite;

         /**
          * @brief Intermediate data for time trigger for state
          */
         int mNextSWrite;
   };

   inline bool TimestepCounter::triggerASCII() const
   {
      return this->mIsATriggered;
   }

   inline bool TimestepCounter::triggerState() const
   {
      return this->mIsSTriggered;
   }

   inline int TimestepCounter::steps() const
   {
      return this->mSteps;
   }

}

#endif // TIMESTEPCOUNTER_HPP
