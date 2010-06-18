/** \file TimestepParameters.hpp
 *  \brief Implementation of a the timestep parameters
 */

#ifndef TIMESTEPPARAMETERS_HPP
#define TIMESTEPPARAMETERS_HPP

// Configuration includes
//
@include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a the timestep parameters
    *
    * \bug Review if all implemented methods are really required (especially the error storage ones)
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
         TimestepParameters(DynamoFloat time, DynamoFloat timestep);

         /**
          * @brief Constructor
          */
         ~TimestepParameters() {};

         /**
          * @brief Has the timestep been computed (important for example predictor/corrector scheme)?
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
         DynamoFloat dt() const;

         /**
          * @brief Get the previous timestep value
          */
         DynamoFloat oldDt() const;

         /**
          * @brief Get the time value
          */
         DynamoFloat time() const;

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
         void updateTimestep(DynamoFloat dt);

         /**
          * @brief Update time
          */
         void updateTime();

         /**
          * @brief Set time
          *
          * @param time Set time
          */
         void setTime(const DynamoFloat time);

         /**
          * @brief Print some information about run
          */
         void printInfo() const;

         /**
          * @brief Update the time derivative truncation error
          *
          * @param order Order of the derivative
          * @param err New value of the error
          */
         void updateDtError(const int order, DynamoFloat err);

         /**
          * @brief Reset the time derivative truncation error
          *
          * @param order Order of the derivative
          */
         void resetDtError(const int order);

         /**
          * @brief Get the time derivative truncation error
          *
          * @param order Order of the derivative
          */
         DynamoFloat dtError(const int order) const;

         /**
          * @brief Get current global error value
          */
         DynamoFloat error() const;

         /**
          * @brief Update global error
          *
          * @param err New error
          */
         void updateError(DynamoFloat err);

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
         DynamoFloat   mTime;

         /**
          * @brief Length of the timestep
          */
         DynamoFloat   mTimestep;

         /**
          * @brief Length of the timestep
          */
         DynamoFloat   mOldTimestep;

         /**
          * @brief Storage for general timestep error
          */
         DynamoFloat mError;

         /**
          * @brief Store for the derivative truncation error
          */
         Array mDtError;

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

   inline DynamoFloat TimestepParameters::time() const
   {
      return this->mTime;
   }

   inline DynamoFloat TimestepParameters::dt() const
   {
      return this->mTimestep;
   }

   inline DynamoFloat TimestepParameters::oldDt() const
   {
      return this->mOldTimestep;
   }

   inline void TimestepParameters::updateTime()
   {
      this->mTime += this->mTimestep;
   }

   inline void TimestepParameters::setTime(const DynamoFloat time)
   {
      this->mTime = time;
   }

   inline DynamoFloat TimestepParameters::error() const
   {
      return this->mError;
   }

   inline DynamoFloat TimestepParameters::dtError(const int order) const
   {
      return this->mDtError(order-1);
   }

   /// Typedef for a shared pointer to a TimestepParameters
   typedef  EPMSHARED_PTR<TimestepParameters>  SmartTStepParameters;

}

#endif // TIMESTEPPARAMETERS_HPP
