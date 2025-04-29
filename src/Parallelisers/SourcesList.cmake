#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      LoadSplitterBase.cpp
      WorkflowBase.cpp
      )
if(EPMDYNAMO_MPI)
   list(APPEND ${SrcSubDir}_Srcs MPIWorkflow.cpp)
else(EPMDYNAMO_MPI)
   list(APPEND ${SrcSubDir}_Srcs SerialWorkflow.cpp)
endif(EPMDYNAMO_MPI)

#
# include source files from subdirectories
#
set(SrcSubSubDirs
      SplittingAlgorithms
      )

foreach(SrcSubSubDir ${SrcSubSubDirs})
   set(SrcSubSubDirPath ${SrcSubDirPath}/${SrcSubSubDir})
   include(${SrcSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_Srcs ${SrcSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubDir)
