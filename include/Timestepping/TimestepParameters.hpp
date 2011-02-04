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
          * @brief Get all the CFLs
          */
         const Array& globalCFLs() const;
         const Array& rtpCFLs() const;
         Array& rRTPCFLs();
         const Array& specCFLs() const;
         Array& rSpecCFLs();
         const Array& errCFLs() const;
         Array& rErrCFLs();
         const Array& rtpCFLPos() const;
         Array& rRTPCFLPos();
         const Array& specCFLPos() const;
         Array& rSpecCFLPos();
         const ArrayI& specCFLDeg() const;
         ArrayI& rSpecCFLDeg();

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
         void setHozRTPVCFL(EPMFloat dt, EPMFloat pos);
         void setHozRTPVBCFL(EPMFloat dt, EPMFloat pos);

         /**
          * @brief Set the local radial RTP CFL
          */
         void setRadRTPVCFL(EPMFloat dt, EPMFloat pos);
         void setRadRTPVBCFL(EPMFloat dt, EPMFloat pos);

         /**
          * @brief Set the local horizontal spectral CFL
          */
         void setTorHozSpecVCFL(EPMFloat dt, EPMFloat pos, int l);
         void setTorHozSpecVBCFL(EPMFloat dt, EPMFloat pos, int l);

         /**
          * @brief Set the local horizontal spectral CFL
          */
         void setPolHozSpecVCFL(EPMFloat dt, EPMFloat pos, int l);
         void setPolHozSpecVBCFL(EPMFloat dt, EPMFloat pos, int l);

         /**
          * @brief Set the local radial spectral CFL
          */
         void setPolRadSpecVCFL(EPMFloat dt, EPMFloat pos, int l);
         void setPolRadSpecVBCFL(EPMFloat dt, EPMFloat pos, int l);

         /**
          * @brief Set the error control CFL
          */
         void setErrorCFL(EPMFloat dt, int cId);

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
          * @brief Storage for the different CFL conditions
          */
         Array mGlobalCFLs;
         Array mRTPCFLs;
         Array mSpecCFLs;
         Array mRTPCFLPos;
         Array mSpecCFLPos;
         ArrayI mSpecCFLDeg;

         /**
          * @brief Storage for the different  error CFL conditions
          */
         Array mErrCFLs;
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
      cfl = std::min(cfl, this->mErrCFLs.minCoeff());

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

   inline const Array& TimestepParameters::specCFLs() const
   {
      return this->mSpecCFLs;
   }

   inline const Array& TimestepParameters::errCFLs() const
   {
      return this->mErrCFLs;
   }

   inline const Array& TimestepParameters::rtpCFLPos() const
   {
      return this->mRTPCFLPos;
   }

   inline const Array& TimestepParameters::specCFLPos() const
   {
      return this->mSpecCFLPos;
   }

   inline const ArrayI& TimestepParameters::specCFLDeg() const
   {
      return this->mSpecCFLDeg;
   }

   inline Array& TimestepParameters::rRTPCFLs()
   {
      return this->mRTPCFLs;
   }

   inline Array& TimestepParameters::rSpecCFLs()
   {
      return this->mSpecCFLs;
   }

   inline Array& TimestepParameters::rErrCFLs()
   {
      return this->mErrCFLs;
   }

   inline Array& TimestepParameters::rRTPCFLPos()
   {
      return this->mRTPCFLPos;
   }

   inline Array& TimestepParameters::rSpecCFLPos()
   {
      return this->mSpecCFLPos;
   }

   inline ArrayI& TimestepParameters::rSpecCFLDeg()
   {
      return this->mSpecCFLDeg;
   }



   inline void TimestepParameters::setInertialCFL(EPMFloat dt)
   {
      this->mGlobalCFLs(0) = dt;
   }

   inline void TimestepParameters::setTorsionalCFL(EPMFloat dt)
   {
      this->mGlobalCFLs(1) = dt;
   }

   inline void TimestepParameters::setHozRTPVCFL(EPMFloat dt, EPMFloat pos)
   {
      this->mRTPCFLs(0) = dt;
      this->mRTPCFLPos(0) = pos;
   }

   inline void TimestepParameters::setHozRTPVBCFL(EPMFloat dt, EPMFloat pos)
   {
      this->mRTPCFLs(1) = dt;
      this->mRTPCFLPos(1) = pos;
   }

   inline void TimestepParameters::setRadRTPVCFL(EPMFloat dt, EPMFloat pos)
   {
      this->mRTPCFLs(2) = dt;
      this->mRTPCFLPos(2) = pos;
   }

   inline void TimestepParameters::setRadRTPVBCFL(EPMFloat dt, EPMFloat pos)
   {
      this->mRTPCFLs(3) = dt;
      this->mRTPCFLPos(3) = pos;
   }

   inline void TimestepParameters::setTorHozSpecVCFL(EPMFloat dt, EPMFloat pos, int l)
   {
      this->mSpecCFLs(0) = dt;
      this->mSpecCFLPos(0) = pos;
      this->mSpecCFLDeg(0) = l;
   }

   inline void TimestepParameters::setTorHozSpecVBCFL(EPMFloat dt, EPMFloat pos, int l)
   {
      this->mSpecCFLs(1) = dt;
      this->mSpecCFLPos(1) = pos;
      this->mSpecCFLDeg(1) = l;
   }

   inline void TimestepParameters::setPolHozSpecVCFL(EPMFloat dt, EPMFloat pos, int l)
   {
      this->mSpecCFLs(2) = dt;
      this->mSpecCFLPos(2) = pos;
      this->mSpecCFLDeg(2) = l;
   }

   inline void TimestepParameters::setPolHozSpecVBCFL(EPMFloat dt, EPMFloat pos, int l)
   {
      this->mSpecCFLs(3) = dt;
      this->mSpecCFLPos(3) = pos;
      this->mSpecCFLDeg(3) = l;
   }

   inline void TimestepParameters::setPolRadSpecVCFL(EPMFloat dt, EPMFloat pos, int l)
   {
      this->mSpecCFLs(4) = dt;
      this->mSpecCFLPos(4) = pos;
      this->mSpecCFLDeg(4) = l;
   }

   inline void TimestepParameters::setPolRadSpecVBCFL(EPMFloat dt, EPMFloat pos, int l)
   {
      this->mSpecCFLs(5) = dt;
      this->mSpecCFLPos(5) = pos;
      this->mSpecCFLDeg(5) = l;
   }

   inline void TimestepParameters::setErrorCFL(EPMFloat dt, int cId)
   {
      this->mErrCFLs(cId) = dt;
   }

   /// Typedef for a shared pointer to a TimestepParameters
   typedef  EPMSHARED_PTR<TimestepParameters>  SmartTStepParameters;

}

#endif // TIMESTEPPARAMETERS_HPP
