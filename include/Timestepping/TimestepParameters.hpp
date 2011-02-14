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

         /**
          * @brief Get the array of global CFLs
          */
         const Array& globalCFLs() const;

         /**
          * @brief Get the array of RTP CFLs
          */
         const Array& rtpCFLs() const;

         /**
          * @brief Get the error controlled CFL
          */
         EPMFloat errCFL() const;

         /**
          * @brief Get the radial position of the RTP CFL condition
          */
         const Array& rtpCFLPos() const;

         /**
          * @brief Set the array of RTP CFL conditions
          */
         Array& rRTPCFLs();

         /**
          * @brief Set the array of RTP CFL conditions positions
          */
         Array& rRTPCFLPos();

         /**
          * @brief Set the inertial global CFL
          */
         void setInertialCFL(EPMFloat dt);

         /**
          * @brief Set the torsional global CFL
          */
         void setTorsionalCFL(EPMFloat dt);

         /**
          * @brief Set the local horizontal RTP CFL
          */
         void setHozRTPCFL(EPMFloat dt, EPMFloat pos);

         /**
          * @brief Set the local radial RTP CFL
          */
         void setRadRTPCFL(EPMFloat dt, EPMFloat pos);

         /**
          * @brief Set the error control CFL
          */
         void setErrorCFL(EPMFloat dt);

         /**
          * @brief Compute the active CFL condition from all the stored values
          */
         EPMFloat getCFL() const;

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

         /**
          * @brief Storage for the global CFL conditions
          */
         Array mGlobalCFLs;

         /**
          * @brief Storage for the RTP CFL conditions
          */
         Array mRTPCFLs;

         /**
          * @brief Storage for the error CFL condition
          */
         EPMFloat mErrCFL;

         /**
          * @brief Storage for the radial positions of the RTP CFL conditions
          */
         Array mRTPCFLPos;
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

   inline EPMFloat TimestepParameters::getCFL() const
   {
      // Get the minimal global CFL
      EPMFloat cfl = this->mGlobalCFLs.minCoeff();

      // Get the minimal local RTP CFL
      cfl = std::min(cfl, this->mRTPCFLs.minCoeff());

      // Get the minimal error control CFL
      cfl = std::min(cfl, this->mErrCFL);

      return cfl;
   }

   inline const Array& TimestepParameters::globalCFLs() const
   {
      return this->mGlobalCFLs;
   }

   inline const Array& TimestepParameters::rtpCFLs() const
   {
      return this->mRTPCFLs;
   }

   inline EPMFloat TimestepParameters::errCFL() const
   {
      return this->mErrCFL;
   }

   inline void TimestepParameters::setInertialCFL(EPMFloat dt)
   {
      if(dt != 0.0)
      {
         this->mGlobalCFLs(0) = dt;
      } else
      {
         this->mGlobalCFLs(0) = 1.0;
      }
   }

   inline void TimestepParameters::setTorsionalCFL(EPMFloat dt)
   {
      if(dt != 0.0)
      {
         this->mGlobalCFLs(1) = dt;
      } else
      {
         this->mGlobalCFLs(1) = 1.0;
      }
   }

   inline void TimestepParameters::setHozRTPCFL(EPMFloat dt, EPMFloat pos)
   {
      this->mRTPCFLs(0) = dt;
      this->mRTPCFLPos(0) = pos;
   }

   inline void TimestepParameters::setRadRTPCFL(EPMFloat dt, EPMFloat pos)
   {
      this->mRTPCFLs(1) = dt;
      this->mRTPCFLPos(1) = pos;
   }

   inline void TimestepParameters::setErrorCFL(EPMFloat dt)
   {
      this->mErrCFL = dt;
   }

   inline const Array& TimestepParameters::rtpCFLPos() const
   {
      return this->mRTPCFLPos;
   }

   inline Array& TimestepParameters::rRTPCFLs()
   {
      return this->mRTPCFLs;
   }

   inline Array& TimestepParameters::rRTPCFLPos()
   {
      return this->mRTPCFLPos;
   }

   /// Typedef for a shared pointer to a TimestepParameters
   typedef  EPMSHARED_PTR<TimestepParameters>  SmartTStepParameters;

}

#endif // TIMESTEPPARAMETERS_HPP
