#
# Create sources list
#
set(${SrcSubDir}_${SrcSubSubDir}_Srcs
      PCTimestepConfig.cpp
      )

#
# include source files from subdirectories
#
set(SrcSubSubSubDirs
      Theta
      )

foreach(SrcSubSubSubDir ${SrcSubSubSubDirs})
   set(SrcSubSubSubDirPath ${SrcSubSubDirPath}/${SrcSubSubSubDir})
   include(${SrcSubSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_${SrcSubSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_${SrcSubSubDir}_Srcs ${SrcSubSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubSubDir)
