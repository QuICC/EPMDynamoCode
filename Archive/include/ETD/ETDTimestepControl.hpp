/** \file ETDTimestepControl.hpp
 */

#ifndef ETDTIMESTEPCONTROL_HPP
#define ETDTIMESTEPCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/Typedefs.hpp"
#include "Equations/EquationParameters.hpp"
#include "Timestepping/TimestepControl.hpp"
#include "Timestepping/TimestepController.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   class ETDTimestepControl: public TimestepControl
   {
      public:
         /**
          * @brief Constructor
          */
         ETDTimestepControl(TimestepParameters &params, const EquationParameters &eqParams);

         /**
          * @brief Destructor
          */
         virtual ~ETDTimestepControl();
         
      protected:

      private:
   };

}

#endif // ETDTIMESTEPCONTROL_HPP
