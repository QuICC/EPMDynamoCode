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
         TimestepController(DynamoFloat beta1, DynamoFloat beta2, DynamoFloat alpha, TimestepParameters& params);
         
         /**
          * @brief Constructor for predefined type
          *
          * @param type Type of the controller
          * @param order Order of the timestep scheme
          * @param params Timestep parameters
          */
         TimestepController(int type, int order, TimestepParameters& params);

         /**
          * @brief Destructor
          */
         virtual ~TimestepController() {};

         /**
          * @brief Elementary controller type
          */
         static const int ELEMENTARY = 1;
 
         /**
          * @brief PI4.2 controller type
          */
         static const int PI42 = 2;        
 
         /**
          * @brief H211b digital filter controller type
          */
         static const int H211B = 3;        

         /**
          * @brief General routine for first order timestep controller (see Soederlind)
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          * @param errn Previous error
          * @param errn_1 Error at step before previous
          */
         DynamoFloat nextTimestep(DynamoFloat errn, DynamoFloat errn_1);

      protected:

      private:
         /**
          * @brief \f$\beta_1\f$ parameter of the elementary controller
          */
         static const DynamoFloat ELEMENTARY_KBETA1;

         /**
          * @brief \f$\beta_2\f$ parameter of the elementary controller
          */
         static const DynamoFloat ELEMENTARY_KBETA2;

         /**
          * @brief \f$\alpha\f$ parameter of the elementary controller
          */
         static const DynamoFloat ELEMENTARY_ALPHA;

         /**
          * @brief \f$\beta_1\f$ parameter of the PI4.2 controller
          */
         static const DynamoFloat PI42_KBETA1;

         /**
          * @brief \f$\beta_2\f$ pf the PI4.2 controller
          */
         static const DynamoFloat PI42_KBETA2;

         /**
          * @brief \f$\alpha\f$ pa the PI4.2 controller
          */
         static const DynamoFloat PI42_ALPHA;

         /**
          * @brief \f$\beta_1\f$ pf the H211b controller
          */
         static const DynamoFloat H211B_KBETA1;

         /**
          * @brief \f$\beta_2\f$ pf the H211b controller
          */
         static const DynamoFloat H211B_KBETA2;

         /**
          * @brief \f$\alpha\f$ pa the H211b controller
          */
         static const DynamoFloat H211B_ALPHA;

         /**
          * @brief First order timestep controller \f$\beta_1\f$ paramter
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          *  - \f$k\beta_1 = 1\f$ for elementary controller
          *  - \f$k\beta_1 = 3/5\f$ for PI4.2 controller
          *  - \f$k\beta_1 = 1/4\f$ for H211b digital filter (b=4)
          */
         DynamoFloat mControllerBeta1;

         /**
          * @brief First order timestep controller \f$\beta_2\f$ paramter
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          *  - \f$k\beta_2 = 0\f$ for elementary controller
          *  - \f$k\beta_2 = -1/5\f$ for PI4.2 controller
          *  - \f$k\beta_2 = 1/4\f$ for H211b digital filter (b=4)
          */
         DynamoFloat mControllerBeta2;

         /**
          * @brief First order timestep controller \f$\alpha\f$ paramter
          *
          * \f[ h_{n+1} = \left(\frac{\epsilon}{\hat{r}_{n}}\right)^{\beta_1}\left(\frac{\epsilon}{\hat{r}_{n-1}}\right)^{\beta_2}\left(\frac{h_{n}}{h_{n-1}}\right)^{-\alpha} h_{n} \f]
          *
          *  - \f$\alpha = 0\f$ for elementary controller
          *  - \f$\alpha = 0\f$ for PI4.2 controller
          *  - \f$\alpha = 1/4\f$ for H211b digital filter (b=4)
          */
         DynamoFloat mControllerAlpha;

         /**
          * @brief Reference to TimestepParameters object
          */
         TimestepParameters& mrParams;

         /**
          * @brief Set the parameters of the timestep controller
          *
          * @param kb1 \f$k\beta_1\f$ parameter
          * @param kb2 \f$k\beta_2\f$ parameter
          * @param a \f$\alpha\f$ parameter
          * @param k Order of the scheme
          */
         void setControllerParameters(DynamoFloat kb1, DynamoFloat kb2, DynamoFloat a, int k);

         /**
          * @brief Set parameters for elementary controller
          *
          * @param k Order of method
          */
         void setElementaryController(int k);

         /**
          * @brief Set parameters for PI4.2 controller
          *
          * @param k Order of method
          */
         void setPI42Controller(int k);

         /**
          * @brief Set parameters for H211b digital filter
          *
          * @param k Order of method
          */
         void setH211bController(int k);
   };

}

#endif // TIMESTEPCONTROLLER_HPP
