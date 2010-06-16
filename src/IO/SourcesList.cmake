#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      FileBase.cpp
      IOSystem.cpp
      )

#
# include source files from subdirectories
#
set(SrcSubSubDirs ASCII HDF5)

foreach(SrcSubSubDir ${SrcSubSubDirs})
   set(SrcSubSubDirPath ${SrcSubDirPath}/${SrcSubSubDir})
   include(${SrcSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_Srcs ${SrcSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubDir)
