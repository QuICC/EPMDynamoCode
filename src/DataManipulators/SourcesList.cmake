#
# Set list of source files
#
set(${SrcSubDir}_Srcs
   SerialFFTManipulator.cpp
      )

if(EPMDYNAMO_MPI)
  list(APPEND ${SrcSubDir}_Srcs  MPIFFTManipulator.cpp)
endif(EPMDYNAMO_MPI)
