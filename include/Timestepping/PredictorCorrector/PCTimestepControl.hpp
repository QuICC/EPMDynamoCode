/** \file PCTimestepControl.hpp
 *  \brief Implementation of a specialised timestep control for predictor/corrector schemes
 */

#ifndef PCTIMESTEPCONTROL_HPP
#define PCTIMESTEPCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Timestepping/TimestepControl.hpp"
#include "Timestepping/TimestepController.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a specialised timestep control for predictor/corrector schemes
    *
    * \bug This needs a big cleaning and restructuring!!!
    */
   class PCTimestepControl: public TimestepControl
   {
      public:
         /**
          * @brief Constructor
          *
          * @param params Timestep parameters
          * @param eqParams Equation parameters
          */
         PCTimestepControl(TimestepParameters &params, const EquationParameters &eqParams);

         /**
          * @brief Constructor
          */
         ~PCTimestepControl() {};

         /**
          * @brief Check that the timestepping is converging
          *
          * @param step ? Current step in a multistep transform ?
          */
         void checkConvergence(int step);

         /**
          * @brief Update the timestep to use
          */
         void updateTimestep();

         /**
          * @brief Reset value of the PC error
          */
         void resetPCError();

      protected:
         /**
          * @brief Counter for the number of Corrector iterations
          */
         int  mCounterCorrIts;

         /**
          * @brief Timestep obtained through the PC error constrains
          */
         EPMFloat   mPCTStep;

         /**
          * @brief Storage for the norm of the current corrector iteration
          */
         EPMFloat   mCorrectionNorm;

         /**
          * @brief Storage for the norm of the previous corrector iteration
          */
         EPMFloat   mOldCorrectionNorm;

         /**
          * @brief Storage for the total error of the PC step
          */
         EPMFloat   mPCError;

         /**
          * @brief Storage for the total error of the previous PC step
          */ 
         EPMFloat   mOldPCError;

         /**
          * @brief Adaptive timestep controller
          */
         TimestepController  mController;

         /**
          * @brief Get PC errors and compute required timestep
          */
         void setErrorConstrains();

         /**
          * @brief Test for PC error timestep
          *
          * @param rDt Timestep lenght old/new
          */
         void testPCTimestep(EPMFloat& rDt);

         /**
          * @brief Compute the timestep resulting from PC errors
          */
         void computePCTimestep();

      private:

         /**
          * @brief Tot amount of correction from corrector
          *
          * \bug this needs to be restructured
          */
         EPMFloat mTotCorr;
   };

}

#endif // PCTIMESTEPCONTROL_HPP
