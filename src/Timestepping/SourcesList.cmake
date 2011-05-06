#
# Set list of source files
#
set(${SrcSubDir}_Srcs
      ErrorControl.cpp
      IterativeSchemeBase.cpp
      RuntimeControl.cpp
      SchemeIteration.cpp
      TimestepConfig.cpp
      TimestepControlBase.cpp
      TimestepController.cpp
      TimestepCounter.cpp
      TimestepParameters.cpp
      TimestepSchemeBase.cpp
      )

#
# include source files from subdirectories
#
set(SrcSubSubDirs
      ErrorNorms
      PredictorCorrector
      Traits
      )

foreach(SrcSubSubDir ${SrcSubSubDirs})
   set(SrcSubSubDirPath ${SrcSubDirPath}/${SrcSubSubDir})
   include(${SrcSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_Srcs ${SrcSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubDir)
