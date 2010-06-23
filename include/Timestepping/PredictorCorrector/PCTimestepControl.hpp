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
#include "Timestepping/TimestepControlBase.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a specialised timestep control for predictor/corrector schemes
    */
   class PCTimestepControl: public TimestepControlBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param tsParams Timestep parameters
          * @param eqParams Equation parameters
          */
         PCTimestepControl(TimestepParameters &tsParams, const EquationParameters &eqParams);

         /**
          * @brief Constructor
          */
         ~PCTimestepControl() {};

         /**
          * @brief Check that the timestepping is converging
          *
          * @param step Current step
          */
         virtual void checkConvergence(int step);

         /**
          * @brief Update the timestep to use
          */
         virtual void updateTimestep();

      protected:
         /**
          * @brief Counter for the number of Corrector iterations
          */
         int  mCorrectorIts;

         /**
          * @brief Storage for the norm of the current corrector iteration
          */
         EPMFloat   mCorrectionNorm;

         /**
          * @brief Storage for the norm of the previous corrector iteration
          */
         EPMFloat   mOldCorrectionNorm;

         /**
          * @brief Get PC errors and compute required timestep
          */
         void setErrorConstrains();

      private:
   };

}

#endif // PCTIMESTEPCONTROL_HPP
