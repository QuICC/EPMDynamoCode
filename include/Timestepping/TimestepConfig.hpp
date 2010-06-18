/** \file TimestepConfig.hpp
 *  \brief Some general timestepping constants
 */

#ifndef TIMESTEPCONFIG_HPP
#define TIMESTEPCONFIG_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Some general timestepping constants independent of the actual scheme
    */
   class TimestepConfig
   {
      public:
         /**
          * @brief Minimum timestep length
          */
         static const DynamoFloat MIN_TIMESTEP;

         /**
          * @brief Maximal value for timestep
          *
          * This is useful in case the adaptative control is not strong enough or 
          * if for any reason a given timestep should be enforced
          */
         static  DynamoFloat MAXIMUM_TIMESTEP;

         /**
          * @brief Courant number
          */
         static const DynamoFloat COURANT_NUMBER;

         /**
          * @brief Adaptive timestep scheme margin \f$p\f$
          *
          * The timestep will only be modified if the new timestep is outside the
          * \f$ [(1.0-p) dt, (1.0+p) dt] \f$
          */
         static const DynamoFloat TIMESTEP_MARGIN;

         /**
          * @brief Possibility to overrule the maximum allowed timestep
          *
          * @param newMax New maximum allowed timestep
          */
         static void overruleMaxTimestep(DynamoFloat newMax);

         /**
          * @brief Multiplicative constant for tolerance ( < 1)
          */
         static const DynamoFloat TIMESTEP_ERROR_THETA;

         /**
          * @brief Maximum error tolerance
          */
         static const DynamoFloat TIMESTEP_ERROR_TOLERANCE;

         /**
          * @brief Maximum error epsilon (theta*tol)
          */
         static const DynamoFloat TIMESTEP_ERROR_EPSILON;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         TimestepConfig();

         /**
          * @brief Destructor
          */
         virtual ~TimestepConfig() {};
   };

}

#endif // TIMESTEPCONFIG_HPP
