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

namespace EPMDynamo {

   TimestepController::TimestepController(DynamoFloat beta1, DynamoFloat beta2, DynamoFloat alpha, TimestepParameters &params)
      : mControllerBeta1(beta1), mControllerBeta2(beta2), mControllerAlpha(alpha), mrParams(params)
   {
   }

   TimestepController::TimestepController(int type, int order, TimestepParameters &params)
      : mControllerBeta1(1.0), mControllerBeta2(0.0), mControllerAlpha(0.0), mrParams(params)
   {
      // Setup an elementary controller of order order
      if(type == ELEMENTARY)
      {
         this->setElementaryController(order);
      // Setup an PI4.2 controller of order order
      } else if(type == PI42)
      {
         this->setPI42Controller(order);
      // Setup a H211b digital filter controller of order order
      } else if(type == H211B)
      {
         this->setH211bController(order);
      }
   }

   DynamoFloat TimestepController::nextTimestep(DynamoFloat errn, DynamoFloat errn_1)
   {
      // Include multiplicative for tolerance
      DynamoFloat epsilon = TimestepConfig::TIMESTEP_ERROR_EPSILON;

      // Get previous timestep
      DynamoFloat hn_1 = this->mrParams.oldDt();
      // Get current timestep
      DynamoFloat hn = this->mrParams.dt();
      // Storage for new timestep
      DynamoFloat hn1 = hn;

      // Compute new step with first order controller
      hn1 *= std::pow(epsilon/errn, this->mControllerBeta1)*std::pow(epsilon/errn_1, this->mControllerBeta2)*std::pow(hn/hn_1, -this->mControllerAlpha);

      return hn1;
   }

   void TimestepController::setControllerParameters(DynamoFloat kb1, DynamoFloat kb2, DynamoFloat a, int k)
   {
      // Set beta_1 parameter
      this->mControllerBeta1 = kb1/static_cast<double>(k);

      // Set beta_2 parameter
      this->mControllerBeta2 = kb2/static_cast<double>(k);

      // Set alpha parameter
      this->mControllerAlpha = a;
   }

   void TimestepController::setElementaryController(int k)
   {
      this->setControllerParameters(ELEMENTARY_KBETA1, ELEMENTARY_KBETA2, ELEMENTARY_ALPHA, k);
   }

   void TimestepController::setPI42Controller(int k)
   {
      this->setControllerParameters(PI42_KBETA1, PI42_KBETA2, PI42_ALPHA, k);
   }

   void TimestepController::setH211bController(int k)
   {
      this->setControllerParameters(H211B_KBETA1, H211B_KBETA2, H211B_ALPHA, k);
   }

   const DynamoFloat TimestepController::ELEMENTARY_KBETA1 = 1.0;

   const DynamoFloat TimestepController::ELEMENTARY_KBETA2 = 0.0;

   const DynamoFloat TimestepController::ELEMENTARY_ALPHA = 0.0;

   const DynamoFloat TimestepController::PI42_KBETA1 = 3.0/5.0;

   const DynamoFloat TimestepController::PI42_KBETA2 = -1.0/5.0;

   const DynamoFloat TimestepController::PI42_ALPHA = 0.0;

   const DynamoFloat TimestepController::H211B_KBETA1 = 1.0/4.0;

   const DynamoFloat TimestepController::H211B_KBETA2 = 1.0/4.0;

   const DynamoFloat TimestepController::H211B_ALPHA = 1.0/4.0;

}
