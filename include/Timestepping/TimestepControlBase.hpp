/** \file TimestepControlBase.hpp
 *  \brief Base of the implementation of a simple timestep control
 */

#ifndef TIMESTEPCONTROLBASE_HPP
#define TIMESTEPCONTROLBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "GeneralFields/RTPField.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/TimestepController.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a simple timestep control
    */
   class TimestepControlBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param tsParams Timestep parameters
          * @param eqParams Equation parameters
          * @param ctrlType Type of the timestep controller
          * @param order Order of the timestep scheme
          */
         TimestepControlBase(TimestepParameters &tsParams, const EquationParameters &eqParams, TimestepCtrlTypes ctrlType, int order);

         /**
          * @brief Constructor
          */
         virtual ~TimestepControlBase() {};

         /**
          * @brief Should the simulation keep running?
          */
         bool keepRunning() const;

         /**
          * @brief Check that the timestepping is converging
          *
          * @param step Current timestep
          */
         virtual void checkConvergence(int step) = 0;

         /**
          * @brief Update the timestep to use
          */
         virtual void updateTimestep() = 0;

         /**
          * @brief Update value of the CFL condition timestep for 
          *          non-magnetic case
          *
          * @param velV Velocity rtp space field
          * @param llFactor Spherical geometry factor
          */
         void updateCFLTimestep(const RTPField &velV, const EPMFloat llFactor);

         /**
          * @brief Update value of the CFL condition timestep
          *
          * @param magB Magnetic rtp space field
          * @param velV Velocity rtp space field
          * @param llFactor Spherical geometry factor
          */
         void updateCFLTimestep(const RTPField &magB, const RTPField &velV, const EPMFloat llFactor);

         /**
          * @brief Print some useful information
          */
         void printInfo() const;

      protected:
         /**
          * @brief Boolean to check if simulation should keep running?
          */
         bool mKeepRunning;

         /**
          * @brief Bool for run timestep run
          */
         bool mNeedInit;

         /**
          * @brief Storage for timestep error
          */
         EPMFloat mError;

         /**
          * @brief Storage for previous timestep error
          */
         EPMFloat mOldError;

         /**
          * @brief Get timestep parameters
          */
         TimestepParameters& rTSParams();

         /**
          * @brief Get equation parameters
          */
         const EquationParameters& eqParams() const;

         /**
          * @brief Reset stored error and update previous error
          */
         void resetError();

         /**
          * @brief Get the Simulation wide CFL condition (MPI communication)
          */
         void getSimulationCFLCondition();

         /**
          * @brief test for initialisation timestep condition
          *
          * @param rDt Timestep length
          */
         void testInitialisation(EPMFloat& rDt);

         /**
          * @brief Test for maximum timestep
          *
          * @param rDt Timestep length
          */
         void testMaximumTimestep(EPMFloat& rDt);

         /**
          * @brief Test for CFL condition timestep
          *
          * @param rDt Timestep length
          */
         void testCFLCondition(EPMFloat& rDt);

         /**
          * @brief Use the adaptive timestep condition
          */
         void useAdaptiveTimestep(EPMFloat& rDt);

         /**
          * @brief Include courant number in obtained timestep
          *
          * @param rDt Timestep length
          */
         void useCourantTimestep(EPMFloat& rDt);

         /**
          * @brief Test for allowed timestep window
          *
          * @param dt Timestep length
          */
         void useWindowedTimestep(EPMFloat dt);

         /**
          * @brief Use new timestep without further tests
          */
         void useTimestep(EPMFloat dt);

      private:
         /**
          * @brief CFL condition imposed timestep
          */
         EPMFloat   mCFLTimestep;

         /**
          * @brief Reference to a TimestepParameters object
          */
         TimestepParameters&  mrTSParams;

         /**
          * @brief Reference to an EquationParameters object
          */
         const EquationParameters&  mrEqParams;

         /**
          * @brief Adaptive timestep controller
          */
         TimestepController   mController;
   };

   inline bool TimestepControlBase::keepRunning() const
   {
      return this->mKeepRunning;
   }

   inline TimestepParameters& TimestepControlBase::rTSParams()
   {
      return this->mrTSParams;
   }

   inline const EquationParameters& TimestepControlBase::eqParams() const
   {
      return this->mrEqParams;
   }

}

#endif // TIMESTEPCONTROLBASE_HPP
