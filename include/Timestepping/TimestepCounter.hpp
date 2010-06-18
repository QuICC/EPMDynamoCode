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
          * @param stateRate HDF5 state file save rate
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
          */
         void increment();

         /**
          * @brief Get the number of steps
          */
         int step() const;

         /**
          * @brief is maxtstep unlimited?
          */
         bool unlimited() const;

         /**
          * @brief Trigger ASCII save?
          */
         bool triggerA() const;

         /**
          * @brief Trigger state save?
          */
         bool triggerS() const;

         /**
          * @brief Print some information about run
          */
         void printInfo() const;
         
      protected:

      private:

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
   };

   inline void TimestepCounter::increment()
   {
      ++this->mSteps;
   }

   inline int TimestepCounter::step() const
   {
      return this->mSteps;
   }

}

#endif // TIMESTEPCOUNTER_HPP
