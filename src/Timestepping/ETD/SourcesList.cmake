#
# Create sources list
#
set(${SrcSubDir}_${SrcSubSubDir}_Srcs
      ETDInfluenceKernel.cpp
      )

#
# include source files from subdirectories
#
set(SrcSubSubSubDirs
      ETD1
      ETD2RK
      Operators
      )

foreach(SrcSubSubSubDir ${SrcSubSubSubDirs})
   set(SrcSubSubSubDirPath ${SrcSubSubDirPath}/${SrcSubSubSubDir})
   include(${SrcSubSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubSubDir}_${SrcSubSubSubDir}_Srcs})
      list(APPEND ${SrcSubSubDir}_Srcs ${SrcSubSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubSubDir)
