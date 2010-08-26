#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      EPMException.cpp
      ExecutionTimer.cpp
      MathConstants.cpp
      TimerBase.cpp
      )
if(EPMDYNAMO_MPI)
  list(APPEND ${SrcSubDir}_Srcs MPITimer.cpp)
else(EPMDYNAMO_MPI)
  list(APPEND ${SrcSubDir}_Srcs SerialTimer.cpp)
endif(EPMDYNAMO_MPI)
