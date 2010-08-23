#
# Add projectors precision test
#
set(PrecTests ${PrecTests}
      prec_AssociatedLegendreProjectors
      prec_WorlandProjectors
      )

#
# Create all the tests
#
foreach(PrecTest ${PrecTests})
   set(SrcsList ${PrecTest} ${All_Srcs} ${Tools_Srcs})
   add_executable(${PrecTest} ${SrcsList})
   add_test(${PrecTest} ${CMAKE_CURRENT_BINARY_DIR}/${PrecTest})
endforeach(PrecTest)
