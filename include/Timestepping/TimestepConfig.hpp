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
         static const EPMFloat MIN_TIMESTEP;

         /**
          * @brief Maximal value for timestep
          *
          * This is useful in case the adaptative control is not strong enough or 
          * if for any reason a given timestep should be enforced
          */
         static  EPMFloat MAXIMUM_TIMESTEP;

         /**
          * @brief Maximum number of initialisation steps
          */
         static const int MAX_INITIALISATION_STEPS = 30;

         /**
          * @brief Courant number
          */
         static const EPMFloat COURANT_NUMBER;

         /**
          * @brief Adaptive timestep scheme margin \f$p\f$
          *
          * This constant can be used to "window out" small changes in the timestep.
          * It creates a window of accepted timestep given by 
          * \f$ [(1.0-p) dt, (1.0+p) dt] \f$
          */
         static const EPMFloat TIMESTEP_MARGIN;

         /**
          * @brief Possibility to overrule the maximum allowed timestep
          *
          * @param newMax New maximum allowed timestep
          */
         static void overruleMaxTimestep(EPMFloat newMax);

         /**
          * @brief Multiplicative constant for tolerance ( < 1)
          *
          * The exact meaning of this values depends on the implemented scheme
          */
         static const EPMFloat TIMESTEP_ERROR_THETA;

         /**
          * @brief Maximum error tolerance
          *
          * The exact meaning of this values depends on the implemented scheme
          */
         static const EPMFloat TIMESTEP_ERROR_TOLERANCE;

         /**
          * @brief Maximum error epsilon (theta*tol)
          *
          * The exact meaning of this values depends on the implemented scheme
          */
         static const EPMFloat TIMESTEP_ERROR_EPSILON;

         /**
          * @brief Scaling factor for mixed absolute-relative error control
          *
          * A value of 1.0 will use an absolute error as long as the norm is << 1. 
          * Set it to 0.0 if only absolute errors should be used. But if this is done, the timestep will be extremely strongly constrained as the fields 
          * grow.
          */
         static const EPMFloat TIMESTEP_RELERROR_SCALING;

         /**
          * @brief Max scaling factor for the timestep
          *
          * the NG_AshCode uses a value of 2.0, while Soederlind's paper uses 1.602.
          *
          */
         static const EPMFloat TIMESTEP_MAX_RATIO;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         TimestepConfig() {};

         /**
          * @brief Destructor
          */
         virtual ~TimestepConfig() {};
   };

}

#endif // TIMESTEPCONFIG_HPP
