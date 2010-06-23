/** \file TimestepParameters.hpp
 *  \brief Implementation of a the timestep parameters
 */

#ifndef TIMESTEPPARAMETERS_HPP
#define TIMESTEPPARAMETERS_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the timestep parameters
    */
   class TimestepParameters
   {
      public:

         /**
          * @brief Constructor
          *
          * @param time Current time
          * @param timestep Timestep length
          */
         TimestepParameters(EPMFloat time, EPMFloat timestep);

         /**
          * @brief Constructor
          */
         ~TimestepParameters() {};

         /**
          * @brief Has the timestep been computed?
          *
          * This flag allows to relatively easily implemented "iterative" timestep schemes like predictor/corrector
          */
         bool isNextStep() const;

         /**
          * @brief Are we using a new timestep length?
          */
         bool hasNewDt() const;

         /**
          * @brief Is the timestep rejected?
          */
         bool isRejected() const;

         /**
          * @brief Get the timestep value
          */
         EPMFloat dt() const;

         /**
          * @brief Get the previous timestep value
          */
         EPMFloat oldDt() const;

         /**
          * @brief Get the time value
          */
         EPMFloat time() const;

         /**
          * @brief Set the timestep rejection status
          */
         void rejectTimestep();

         /**
          * @brief Unset the timestep rejection status
          */
         void acceptTimestep();

         /**
          * @brief Toggle the status of the timestep
          */
         void toggleTimestepStatus();

         /**
          * @brief Update timestep length
          *
          * @param dt New timestep length
          */
         void updateTimestep(EPMFloat dt);

         /**
          * @brief Update time
          */
         void updateTime();

         /**
          * @brief Set time
          *
          * @param time Set time
          */
         void setTime(const EPMFloat time);

         /**
          * @brief Print some information about run
          */
         void printInfo() const;

         /**
          * @brief Get current global error value
          */
         EPMFloat error() const;

         /**
          * @brief Update global error
          *
          * @param err New error
          */
         void updateError(EPMFloat err);

         /**
          * @brief Reset error to zero
          */
         void resetError();

      protected:

      private:

         /**
          * @brief Can compute next timestep?
          */
         bool mIsNextStep;

         /**
          * @brief Did timestep dt change ?
          */
         bool mHasNewDt;

         /**
          * @brief Is timestep rejected ?
          */
         bool mIsRejected;

         /**
          * @brief Simulation time
          */
         EPMFloat   mTime;

         /**
          * @brief Length of the timestep
          */
         EPMFloat   mTimestep;

         /**
          * @brief Length of the timestep
          */
         EPMFloat   mOldTimestep;

         /**
          * @brief Storage for general timestep error
          */
         EPMFloat mError;

         /**
          * @brief Check initialisation values
          */
         void checkValues();
   };

   inline bool TimestepParameters::isNextStep() const
   {
      return this->mIsNextStep;
   }

   inline void TimestepParameters::toggleTimestepStatus()
   {
      this->mIsNextStep = !this->mIsNextStep;
   }

   inline bool TimestepParameters::hasNewDt() const
   {
      return this->mHasNewDt;
   }

   inline bool TimestepParameters::isRejected() const
   {
      return this->mIsRejected;
   }

   inline void TimestepParameters::rejectTimestep()
   {
      this->mIsRejected = true;
   }

   inline void TimestepParameters::acceptTimestep()
   {
      this->mIsRejected = false;
   }

   inline EPMFloat TimestepParameters::time() const
   {
      return this->mTime;
   }

   inline EPMFloat TimestepParameters::dt() const
   {
      return this->mTimestep;
   }

   inline EPMFloat TimestepParameters::oldDt() const
   {
      return this->mOldTimestep;
   }

   inline void TimestepParameters::updateTime()
   {
      this->mTime += this->mTimestep;
   }

   inline void TimestepParameters::setTime(const EPMFloat time)
   {
      this->mTime = time;
   }

   inline EPMFloat TimestepParameters::error() const
   {
      return this->mError;
   }

   /// Typedef for a shared pointer to a TimestepParameters
   typedef  EPMSHARED_PTR<TimestepParameters>  SmartTStepParameters;

}

#endif // TIMESTEPPARAMETERS_HPP
