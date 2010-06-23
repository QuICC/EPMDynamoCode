/** \file SimulationControlBase.hpp
 *  \brief Traits independent part of the simulation control
 */

#ifndef SIMULATIONCONTROLBASE_HPP
#define SIMULATIONCONTROLBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/TimestepCounter.hpp"
#include "Timestepping/RuntimeControl.hpp"
#include "IO/ASCII/ControlInterface.hpp"

namespace EPMDynamo {

   /**
    * @brief Base class for the SimulationControl to simplify template implementation
    */
   class SimulationControlBase
   {
      public:
         /**
          * @brief Constructor
          *
          * \param t Initial time
          * \param dt Reached timestep
          * \param maxtstep Maximum number of timesteps
          * \param aRate ASCII save rate
          * \param sRate State file save rate
          * \param wall Wall time
          */
         SimulationControlBase(EPMFloat t, EPMFloat dt, int maxtstep, int aRate, int sRate, EPMFloat wall);

         /**
          * @brief Constructor 
          *
          * \param time Time related values: (0) Initial time, (1) Reached timestep
          * \param runI Run related values (integer type): (0) Maximum number of timsteps, (1) ASCII save rate, (2) State file save rate
          * \param run Run related values (float type): (0) Wall time
          */
         SimulationControlBase(const Array& time, const ArrayI& runI, const Array& run);

         /**
          * @brief Destructor
          */
         virtual ~SimulationControlBase() {};

         /**
          * @brief Last minute initialisation operations before starting time integration
          */
         void preRun();

         /**
          * @brief Early finalisation operations just after ending time integration
          */
         void postRun();

         /**
          * @brief Print post run information
          */
         void printInfo() const;

         /**
          * @brief Get the time integration parameters
          */
         const TimestepParameters&  tsParams() const;

         /**
          * @brief Get the time integration parameters
          */
         TimestepParameters&  tsParams();

         /**
          * @brief Get the time integration counters
          */
         const TimestepCounter&  tsCounter() const;

         /**
          * @brief Get the time integration counters
          */
         TimestepCounter&  tsCounter();

         /**
          * @brief Get the time integration run control
          */
         const RuntimeControl&  runControl() const;

         /**
          * @brief Get the time integration run control
          */
         RuntimeControl&  runControl();

         /**
          * @brief Get the external control interface
          */
         const ControlInterface&  ctrlIface() const;

         /**
          * @brief Get the external control interface
          */
         ControlInterface&  ctrlIface();

         /**
          * @brief Should the simulation keep running?
          */
         bool keepRunning() const;

         /**
          * @brief Check for command from external control interface input
          */
         void checkInterface();
         
      protected:

         /**
          * @brief keep it running?
          */
         bool mKeepRunning;

         /**
          * @brief Timestep parameters object
          */
         TimestepParameters   mTSParams;

         /**
          * @brief Step counter
          */
         TimestepCounter mTSCounter;

         /**
          * @brief Runtime control
          */
         RuntimeControl mRunControl;

         /**
          * @brief Bool to choose between wall time limit and maximum number of timesteps
          */
         const bool  mUseWall;

         /**
          * @brief External control interface
          */
         ControlInterface  mCtrlIface;

      private:
   };

   inline bool SimulationControlBase::keepRunning() const
   {
      return this->mKeepRunning;
   }

   inline const TimestepParameters& SimulationControlBase::tsParams() const
   {
      return this->mTSParams;
   }

   inline TimestepParameters& SimulationControlBase::tsParams()
   {
      return this->mTSParams;
   }

   inline const TimestepCounter& SimulationControlBase::tsCounter() const
   {
      return this->mTSCounter;
   }

   inline TimestepCounter& SimulationControlBase::tsCounter()
   {
      return this->mTSCounter;
   }

   inline const RuntimeControl& SimulationControlBase::runControl() const
   {
      return this->mRunControl;
   }

   inline RuntimeControl& SimulationControlBase::runControl()
   {
      return this->mRunControl;
   }

   inline const ControlInterface& SimulationControlBase::ctrlIface() const
   {
      return this->mCtrlIface;
   }

   inline ControlInterface& SimulationControlBase::ctrlIface()
   {
      return this->mCtrlIface;
   }
}

#endif // SIMULATIONCONTROLBASE_HPP
