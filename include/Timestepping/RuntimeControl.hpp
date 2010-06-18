/** \file RuntimeControl.hpp
 *  \brief Implementation of a simple run time control
 */

#ifndef RUNTIMECONTROL_HPP
#define RUNTIMECONTROL_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

#ifdef EPMDYNAMO_MPI
   #include "General/MPITimer.hpp"
#else
   #include "General/SerialTimer.hpp"
#endif // EPMDYNAMO_MPI

namespace EPMDynamo {

   /**
    * @brief Implementation a small controler dealing with the run time information
    */
   class RuntimeControl: public TimerType
   {
      public:
         /**
          * @brief Constructor
          *
          * @param wall Wall time
          */
         RuntimeControl(EPMFloat wall);

         /**
          * @brief Destructor
          */
         virtual ~RuntimeControl() {};

         /**
          * @brief Check wall time
          */
         bool keepRunning() const;

         /**
          * @brief Get current runtime
          */
         EPMFloat runtime() const;

         /**
          * @brief Update run time
          */
         void  update();

         /**
          * @brief Print some run info
          *
          * @param steps Number of steps performed
          */
         void printInfo(const int steps) const;
         
      protected:

      private:
         /**
          * @brief Wall time
          */
         const EPMFloat mWallTime;

         /**
          * @brief Run time
          */
         EPMFloat mRuntime;
   };

   inline EPMFloat RuntimeControl::runtime() const
   {
      return this->mRuntime;
   }

}

#endif // RUNTIMECONTROL_HPP
