/** \file MPIWorkflow.cpp
 *  \brief Source of the implementation of an MPI workflow control
 */

// System includes
//

// External includes
//
#include <mpi.h>

// Class include
//
#include "Parallelisers/MPIWorkflow.hpp"

// Project includes
//

namespace EPMDynamo {

   MPIWorkflow::MPIWorkflow(const int nCore)
      : WorkflowBase(nCore)
   {
      this->setup();
   }

   void MPIWorkflow::init()
   {
      // Initiliase
      MPI_Init(0, 0);
   }

   void MPIWorkflow::setup()
   {
      // Get MPI size
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      this->checkCores(size);

      // Get MPI rank
      MPI_Comm_rank(MPI_COMM_WORLD, &this->mCoreID);
   }

   void MPIWorkflow::synchronize()
   {
      MPI_Barrier(MPI_COMM_WORLD);
   }

   void MPIWorkflow::finalize()
   {
      // Make sure all finished
      MPI_Barrier(MPI_COMM_WORLD);

      // Finalize MPI system
      MPI_Finalize();
   }

   int MPIWorkflow::rank()
   {
      int rank;

      // Get MPI rank
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      return rank;
   }

}
