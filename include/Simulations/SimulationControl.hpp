/** \file SimulationControl.hpp
 *  \brief Implementation of a general simulation control structure
 */

#ifndef SIMULATIONCONTROL_HPP
#define SIMULATIONCONTROL_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/Parameters/EquationParameters.hpp"
#include "Simulations/SimulationControlBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of simulation controller
    */
   class SimulationControl: public SimulationControlBase
   {
      public:
         /// Typedef for the Timestep control type
         typedef  SimulationConfig::TimestepTraits::TimestepControl   TimestepControl;

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
          * \param maxtime Max integration time
          */
         SimulationControl(EPMFloat t, EPMFloat dt, const EquationParameters& eqParams, int maxtstep, int aRate, int sRate, EPMFloat wall, EPMFloat maxtime);

         /**
          * @brief Constructor
          *
          * \param time Time related values: (0) Initial time, (1) Reached timestep
          * \param eqParams Equation parameters
          * \param runI Run related values (integer type): (0) Maximum number of timsteps, (1) ASCII save rate, (2) State file save rate
          * \param run Run related values (float type): (0) Wall time, (1) Max integration time
          */
         SimulationControl(const Array& time, const EquationParameters& eqParams, const ArrayI& runI, const Array& run);

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

   inline SimulationControl::TimestepControl& SimulationControl::tsControl()
   {
      return this->mTSControl;
   }
}

#endif // SIMULATIONCONTROL_HPP
