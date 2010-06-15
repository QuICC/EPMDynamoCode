/** \file Parallelisation.h
 *  \brief Macros/includes depending on parallelisation
 */

#ifndef PARALLELISATION_H
#define PARALLELISATION_H

#ifdef EPMDYNAMO_MPI
   #include <mpi.h>
   #include "Parallelisers/MPIWorkflow.hpp"
   /// MPI initialisation operations
   #define  EPMDYNAMO_INITIALISER   MPIWorkflow::init()
   /// MPI finalisation operations
   #define  EPMDYNAMO_FINALIZER   MPIWorkflow::finalize()
#else
   #include "Parallelisers/SerialWorkflow.hpp"
   /// Serial initialisation operations
   #define  EPMDYNAMO_INITIALISER   SerialWorkflow::init()
   /// Serial finalisation operations
   #define  EPMDYNAMO_FINALIZER   SerialWorkflow::finalize()
#endif // EPMDYNAMO_MPI

#endif // PARALLELISATION_H
