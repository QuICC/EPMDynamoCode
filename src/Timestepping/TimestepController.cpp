/** \file TimestepController.cpp
 *  \brief Source of the implementation of a set of timestep controllers
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/TimestepController.hpp"

// Project includes
//
#include "Timestepping/TimestepConfig.hpp"
#include "Timestepping/Traits/TimestepControllerTraits.hpp"

namespace EPMDynamo {

   TimestepController::TimestepController(EPMFloat beta1, EPMFloat beta2, EPMFloat alpha, TimestepParameters &params)
      : mControllerBeta1(beta1), mControllerBeta2(beta2), mControllerAlpha(alpha), mrParams(params)
   {
   }

   TimestepController::TimestepController(TimestepCtlTypes type, int order, TimestepParameters &params)
      : mControllerBeta1(1.0), mControllerBeta2(0.0), mControllerAlpha(0.0), mrParams(params)
   {
      // Setup a predefined controller of order order
      this->setPredefinedController(type, order);
   }

   EPMFloat TimestepController::nextTimestep(EPMFloat errn, EPMFloat errn_1)
   {
      // Include multiplicative for tolerance
      EPMFloat epsilon = TimestepConfig::TIMESTEP_ERROR_EPSILON;

      // Get previous timestep
      EPMFloat hn_1 = this->mrParams.oldDt();
      // Get current timestep
      EPMFloat hn = this->mrParams.dt();
      // Storage for new timestep
      EPMFloat hn1 = hn;

      // Compute new step with first order controller
      hn1 *= std::pow(epsilon/errn, this->beta1())*std::pow(epsilon/errn_1, this->beta2())*std::pow(hn/hn_1, -this->alpha());

      return hn1;
   }

   void TimestepController::setControllerParameters(EPMFloat kb1, EPMFloat kb2, EPMFloat a, int k)
   {
      // Set beta_1 parameter
      this->mControllerBeta1 = kb1/static_cast<double>(k);

      // Set beta_2 parameter
      this->mControllerBeta2 = kb2/static_cast<double>(k);

      // Set alpha parameter
      this->mControllerAlpha = a;
   }

   void TimestepController::setPredefinedController(TimestepCtlTypes type, int k)
   {
      // Setup elementary controller
      if(type == ElementaryCtl)
      {
         this->setControllerParameters(TimestepControllerTraits<ElementaryCtl>::KBETA1, TimestepControllerTraits<ElementaryCtl>::KBETA2, TimestepControllerTraits<ElementaryCtl>::ALPHA, k);

      // Setup PI42 controller
      } else if(type == PI42Ctl)
      {
         this->setControllerParameters(TimestepControllerTraits<PI42Ctl>::KBETA1, TimestepControllerTraits<PI42Ctl>::KBETA2, TimestepControllerTraits<PI42Ctl>::ALPHA, k);

      // Setup H211B controller
      } else if(type == H211BCtl)
      {
         this->setControllerParameters(TimestepControllerTraits<H211BCtl>::KBETA1, TimestepControllerTraits<H211BCtl>::KBETA2, TimestepControllerTraits<H211BCtl>::ALPHA, k);

      // Should never happen
      } else
      {
         assert(false);
      }
   }
}
