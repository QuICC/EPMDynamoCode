#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      SimulationBase.cpp
      SimulationConstants.cpp
      SimulationControl.cpp
      SimulationControlBase.cpp
      )

#
# include source files from subdirectories
#
set(SrcSubSubDirs
      Implementations
      )

foreach(SrcSubSubDir ${SrcSubSubDirs})
   set(SrcSubSubDirPath ${SrcSubDirPath}/${SrcSubSubDir})
   include(${SrcSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_Srcs ${SrcSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubDir)
