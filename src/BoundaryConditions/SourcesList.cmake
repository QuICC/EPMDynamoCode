#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      BoundaryCondition.cpp
      HomogeneousBC.cpp
      L0HarmonicBC.cpp
      TimeDependentBC.cpp
      )

#
# include source files from subdirectories
#
set(SrcSubSubDirs
      Homogeneous
      Inhomogeneous
      L0Harmonic
      TimeDependent
      )

foreach(SrcSubSubDir ${SrcSubSubDirs})
   set(SrcSubSubDirPath ${SrcSubDirPath}/${SrcSubSubDir})
   include(${SrcSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_Srcs ${SrcSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubDir)
