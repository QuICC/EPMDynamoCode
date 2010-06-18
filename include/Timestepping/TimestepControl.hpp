/** \file TimestepControl.hpp
 *  \brief Implementation of a simple timestep control
 */

#ifndef TIMESTEPCONTROL_HPP
#define TIMESTEPCONTROL_HPP

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

namespace EPMDynamo {

   /**
    * \brief Implementation of a simple timestep control
    *
    * \bug Review the different timestep tests
    * \bug Improve MPI communication
    */
   class TimestepControl
   {
      public:
         /**
          * @brief Constructor
          *
          * @param params Timestep parameters
          * @param eqParams Equation parameters
          */
         TimestepControl(TimestepParameters &params, const EquationParameters &eqParams);

         /**
          * @brief Constructor
          */
         ~TimestepControl() {};

         /**
          * @brief Should the simulation keep running?
          */
         bool keepRunning() const;

         /**
          * @brief Check that the timestepping is converging
          */
         void checkConvergence(int step);

         /**
          * @brief Update the timestep to use
          */
         void updateTimestep();

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
          * @brief Reference to a TimestepParameters object
          */
         TimestepParameters&  mrParams;

         /**
          * @brief Reference to an EquationParameters object
          */
         const EquationParameters&  mrEqParams;

         /**
          * @brief CFL condition imposed timestep
          */
         EPMFloat   mCFLTimestep;

         /**
          * @brief Get the Simulation wide CFL condition (MPI communication)
          */
         void getSimulationCFLCondition();

         /**
          * @brief test for initialisation timestep
          *
          * @param rDt Timestep length
          */
         void testInitialisationTimestep(EPMFloat& rDt);

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
         void testCFLTimestep(EPMFloat& rDt);

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
         void setWindowedTimestep(EPMFloat dt);

      private:
   };

   inline bool TimestepControl::keepRunning() const
   {
      return this->mKeepRunning;
   }

}

#endif // TIMESTEPCONTROL_HPP
