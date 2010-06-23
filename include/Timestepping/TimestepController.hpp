/** \file TimestepController.hpp
 *  \brief Implementation of a set of timestep controllers
 */

#ifndef TIMESTEPCONTROLLER_HPP
#define TIMESTEPCONTROLLER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/Traits/TimestepControllerTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a set of timestep controllers of different complexities
    */
   class TimestepController
   {
      public:
         /**
          * @brief Constructor with complete flexibility
          *
          * @param beta1 Beta1 parameter of controller
          * @param beta2 Beta2 parameter of controller
          * @param alpha Alpha parameter of controller
          * @param params Timestep parameters
          */
         TimestepController(EPMFloat beta1, EPMFloat beta2, EPMFloat alpha, TimestepParameters& params);
         
         /**
          * @brief Constructor for predefined type
          *
          * @param type Type of the controller
          * @param order Order of the timestep scheme
          * @param params Timestep parameters
          */
         TimestepController(TimestepCtlTypes type, int order, TimestepParameters& params);

         /**
          * @brief Destructor
          */
         virtual ~TimestepController() {};

         /**
          * @brief General routine for first order timestep controller (see Soederlind)
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          * @param errn Previous error
          * @param errn_1 Error at step before previous
          */
         EPMFloat nextTimestep(EPMFloat errn, EPMFloat errn_1);

      protected:

      private:
         /**
          * @brief First order timestep controller \f$\beta_1\f$ paramter
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          *  - \f$k\beta_1 = 1\f$ for elementary controller
          *  - \f$k\beta_1 = 3/5\f$ for PI4.2 controller
          *  - \f$k\beta_1 = 1/4\f$ for H211b digital filter (b=4)
          */
         EPMFloat mControllerBeta1;

         /**
          * @brief First order timestep controller \f$\beta_2\f$ paramter
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          *  - \f$k\beta_2 = 0\f$ for elementary controller
          *  - \f$k\beta_2 = -1/5\f$ for PI4.2 controller
          *  - \f$k\beta_2 = 1/4\f$ for H211b digital filter (b=4)
          */
         EPMFloat mControllerBeta2;

         /**
          * @brief First order timestep controller \f$\alpha\f$ paramter
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          *  - \f$\alpha = 0\f$ for elementary controller
          *  - \f$\alpha = 0\f$ for PI4.2 controller
          *  - \f$\alpha = 1/4\f$ for H211b digital filter (b=4)
          */
         EPMFloat mControllerAlpha;

         /**
          * @brief Reference to TimestepParameters object
          */
         TimestepParameters& mrParams;

         /**
          * @brief Get \f$\beta_1\f$ paramter
          */
         EPMFloat beta1() const;

         /**
          * @brief Get \f$\beta_2\f$ paramter
          */
         EPMFloat beta2() const;

         /**
          * @brief Get \f$\alpha\f$ paramter
          */
         EPMFloat alpha() const;

         /**
          * @brief Set the parameters of the timestep controller
          *
          * @param kb1 \f$k\beta_1\f$ parameter
          * @param kb2 \f$k\beta_2\f$ parameter
          * @param a \f$\alpha\f$ parameter
          * @param k Order of the scheme
          */
         void setControllerParameters(EPMFloat kb1, EPMFloat kb2, EPMFloat a, int k);

         /**
          * @brief Set the parameters for a predefined controller type
          *
          * @param type Type of the controller
          * @param k Order of the controller
          */
          void setPredefinedController(TimestepCtlTypes type, int k);
   };

   inline EPMFloat   TimestepController::beta1() const
   {
      return this->mControllerBeta1;
   }

   inline EPMFloat   TimestepController::beta2() const
   {
      return this->mControllerBeta2;
   }

   inline EPMFloat   TimestepController::alpha() const
   {
      return this->mControllerAlpha;
   }

}

#endif // TIMESTEPCONTROLLER_HPP
