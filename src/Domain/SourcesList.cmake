#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      CoreTruncationBase.cpp
      CoreTruncation.cpp 
      LocalTruncation.cpp
      SimulationTruncation.cpp
      Truncation.cpp
      )

#
# include source files from subdirectories
#
set(SrcSubSubDirs
      Core
      Simulation
      )

foreach(SrcSubSubDir ${SrcSubSubDirs})
   set(SrcSubSubDirPath ${SrcSubDirPath}/${SrcSubSubDir})
   include(${SrcSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_Srcs ${SrcSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubDir)
