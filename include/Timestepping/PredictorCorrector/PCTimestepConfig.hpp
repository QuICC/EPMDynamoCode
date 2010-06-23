/** \file PCTimestepConfig.hpp
 *  \brief Implementation of specialised timestep config for the predictor/corrector schemes
 */

#ifndef PCTIMESTEPCONFIG_HPP
#define PCTIMESTEPCONFIG_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of specialised timestep control for the predictor/corrector schemes
    */
   class PCTimestepConfig
   {
      public:
         /**
          * @brief Maximum number of corrector iterations
          */
         static const int MAX_CORRECTOR_ITERATIONS = 10;

         /**
          * @brief Maximum number of initialisation steps
          */
         static const int MAX_INITIALISATION_STEPS = 30;

         /**
          * @brief Maximum corrector norm value
          */
         static const EPMFloat MAX_CORRECTION_NORM;

         /**
          * @brief Maximum corrector error value
          */
         static const EPMFloat MAX_CORRECTOR_ERROR;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         PCTimestepConfig() {};

         /**
          * @brief Destructor
          */
         virtual ~PCTimestepConfig() {};
   };

}

#endif // PCTIMESTEPCONFIG_HPP
