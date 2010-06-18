/** \file MPIWorkflow.hpp
 *  \brief Implementation of an MPI workflow control
 */

#ifndef MPIWORKFLOW_HPP
#define MPIWORKFLOW_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Parallelisers/WorkflowBase.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the workflow for an MPI Parallel code
    */
   class MPIWorkflow: public WorkflowBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nCore Number of cores
          */
         MPIWorkflow(const int nCore);

         /**
          * @brief Destructor
          */
         virtual ~MPIWorkflow() {};

         /**
          * @brief Initialise the MPI system
          */
         static void init();

         /**
          * @brief Get the mpi rank
          *
          * This routine should only be used while the Truncation object is not available...
          */
         static int rank();

         /**
          * @brief Finalize the MPI system
          */
         static void finalize();
         
      protected:
         /**
          * @brief Setup the MPI system
          */
         void setup();

      private:
   };

}

#endif // MPIWORKFLOW_HPP
