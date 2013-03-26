#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      )

if(EPMDYNAMO_MPI)
  list(APPEND ${SrcSubDir}_Srcs  MPIFFTManipulator.cpp)
else(EPMDYNAMO_MPI)
  list(APPEND ${SrcSubDir}_Srcs SerialFFTManipulator.cpp)
endif(EPMDYNAMO_MPI)
