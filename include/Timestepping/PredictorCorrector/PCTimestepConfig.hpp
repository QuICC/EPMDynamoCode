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
    *
    * \bug Needs to be reviewed and cleaned
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
          * @brief Maximum timestep autorised by stability
          */
         static const DynamoFloat MAX_CORRECTOR_TSTEP;

         /**
          * @brief Maximum corrector norm value
          */
         static const DynamoFloat MAX_CORRECTION_NORM;

         /**
          * @brief Use the last corrector correction as error diagnostic
          */
         static const int CORRECTOR_NORM = 1;

         /**
          * @brief Use the sum of the corrector correction as error diagnostic
          *
          * This is most of the time the same as the previous since there
          * is normally only on corrector iteration
          */
         static const int SUMMED_CORRECTOR_NORM = 2;

         /**
          * @brief Set the type of error diagnostic to use
          */
         static const int PC_ERROR_TYPE = PCTimestepConfig::SUMMED_CORRECTOR_NORM;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         PCTimestepConfig();

         /**
          * @brief Destructor
          */
         virtual ~PCTimestepConfig() {};
   };

}

#endif // PCTIMESTEPCONFIG_HPP
