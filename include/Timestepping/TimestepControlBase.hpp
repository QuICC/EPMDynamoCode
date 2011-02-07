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
#include "GeneralScalars/SpectralSHScalar.hpp"
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
          * @brief Update value of the RTP CFL conditions timestep for 
          *          non-magnetic case
          *
          * @param velV Velocity rtp space field
          */
         void updateRTPCFLTimestep(const RTPField &velV);

         /**
          * @brief Update value of the RTP CFL condition timestep
          *
          * @param magB Magnetic rtp space field
          * @param velV Velocity rtp space field
          */
         void updateRTPCFLTimestep(const RTPField &magB, const RTPField &velV);

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
         EPMFloat getSimulationCFLCondition();

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
          * @brief Use the adaptive timestep error controller condition
          */
         void useErrorCtrlTimestep();

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

         TimestepController   mCtrlPI42;
         TimestepController   mCtrlH211B;
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
