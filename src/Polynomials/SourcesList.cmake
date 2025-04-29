#
# Set list of general source files
#
set(${SrcSubDir}_Srcs
      AssocLegendreBasis.cpp
      AssociatedLegendrePolynomial.cpp
      AssociatedLegendreOperatorBase.cpp
      AssociatedLegendreOperator.cpp
      ChebyshevPolynomial.cpp
      OnesidedJacobi.cpp
      Polynomial.cpp
      QuadratureRule.cpp
      WeightedPolynomial.cpp
      WorlandPolynomial.cpp
      WorlandChebyshevTraits.cpp
      WorlandLegendreTraits.cpp
      )
#
# include source files from subdirectories
#
set(SrcSubSubDirs
      PolyOps
      )

foreach(SrcSubSubDir ${SrcSubSubDirs})
   set(SrcSubSubDirPath ${SrcSubDirPath}/${SrcSubSubDir})
   include(${SrcSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_Srcs ${SrcSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubDir)
