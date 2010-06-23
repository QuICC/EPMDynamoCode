/** \file SerialWorkflow.hpp
 *  \brief Implementation of the serial workflow control
 */

#ifndef SERIALWORKFLOW_HPP
#define SERIALWORKFLOW_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Parallelisers/WorkflowBase.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the workflow for a serial code
    */
   class SerialWorkflow: public WorkflowBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nCore Number of cores
          */
         SerialWorkflow(const int nCore);

         /**
          * @brief Destructor
          */
         virtual ~SerialWorkflow() {};

         /**
          * @brief Initialise the Serial system
          */
         static void init();

         /**
          * @brief Synchronize
          */
         static void synchronize();
  
         /**
          * @brief Finalize the Serial system
          */
         static void finalize();
         
      protected:
         /**
          * @brief Setup the Serial system
          */
         void setup();

      private:
   };

}

#endif // SERIALWORKFLOW_HPP
