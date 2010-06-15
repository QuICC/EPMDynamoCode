/** \file SimulationControl.hpp
 *  \brief Implementation of a general simulation control structure
 */

#ifndef SIMULATIONCONTROL_HPP
#define SIMULATIONCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Simulations/SimulationControlBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of simulation controller
    *
    * \tparam TSimType Type of simulation
    */
   template <typename TSimType> class SimulationControl: public SimulationControlBase
   {
      public:
         /// Typedef for the Timestep control type
         typedef  typename SimulationTraits<TSimType>::TimestepTraits::TimestepControl   TimestepControl;

         /// Typedef for the equation parameters type
         typedef  typename TSimType::EqParamsType   EqParamsType;

         /**
          * @brief Constructor
          *
          * \param t Initial time
          * \param dt Reached timestep
          * \param eqParams Equation parameters
          * \param maxtstep Maximum number of timesteps
          * \param aRate ASCII save rate
          * \param sRate State file save rate
          * \param wall Wall time
          */
         SimulationControl(DynamoFloat t, DynamoFloat dt, const EqParamsType& eqParams, int maxtstep, int aRate, int sRate, DynamoFloat wall);

         /**
          * @brief Constructor
          *
          * \param time Time related values: (0) Initial time, (1) Reached timestep
          * \param eqParams Equation parameters
          * \param runI Run related values (integer type): (0) Maximum number of timsteps, (1) ASCII save rate, (2) State file save rate
          * \param run Run related values (float type): (0) Wall time
          */
         SimulationControl(const Array& time, const EqParamsType& eqParams, const ArrayI& runI, const Array& run);

         /**
          * @brief Destructor
          */
         virtual ~SimulationControl() {};

         /**
          * @brief Get the TimestepControl object
          */
         TimestepControl&  tsControl();

         /**
          * @brief Update status
          */
         void update();
         
      protected:

         /**
          * @brief Timestep control object
          */
         TimestepControl   mTSControl;

      private:
   };

   template <typename TSimType> inline typename SimulationControl<TSimType>::TimestepControl& SimulationControl<TSimType>::tsControl()
   {
      return this->mTSControl;
   }

   template <typename TSimType> SimulationControl<TSimType>::SimulationControl(DynamoFloat t, DynamoFloat dt, const EqParamsType& eqParams, int maxtstep, int aRate, int sRate, DynamoFloat wall)
      : SimulationControlBase(t, dt, maxtstep, aRate, sRate, wall), mTSControl(mTSParams, eqParams)
   {
   }

   template <typename TSimType> SimulationControl<TSimType>::SimulationControl(const Array& time, const EqParamsType& eqParams, const ArrayI& runI, const Array& run)
      : SimulationControlBase(time, runI, run), mTSControl(mTSParams, eqParams)
   {
   }

   template <typename TSimType> void SimulationControl<TSimType>::update()
   {
      // Check convergence of time integration
      this->tsControl().checkConvergence(this->tsCounter().step());

      if(this->tsParams().isNextStep())
      {
         // Update simulation time
         this->tsParams().updateTime();

         // increment the number of performed timesteps
         this->tsCounter().increment();

         // Update runtime
         this->runControl().update();
      }

      // Keep running simulation?
      if(this->mUseWall)
      {
         this->mKeepRunning = this->tsControl().keepRunning() && this->runControl().keepRunning() && this->ctrlIface().keepRunning();
      } else
      {
         this->mKeepRunning = this->tsControl().keepRunning() && this->tsCounter().keepRunning() && this->ctrlIface().keepRunning();
      }
   }
}

#endif // SIMULATIONCONTROL_HPP
