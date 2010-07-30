#
# Add transform precision test
#
set(PrecTests ${PrecTests}
      prec_AssocLegendreBasis
      prec_FFTransform
      prec_LegendreTransform
      prec_RTPOperators
      prec_SpectralSHTransform
      prec_WorlandBasis
      prec_WorlandTransform
      )

#
# Create all the tests
#
foreach(PrecTest ${PrecTests})
   set(SrcsList ${PrecTest} ${All_Srcs} ${Tools_Srcs})
   add_executable(${PrecTest} ${SrcsList})
   add_test(${PrecTest} ${CMAKE_CURRENT_BINARY_DIR}/${PrecTest})
endforeach(PrecTest)

#
# Create test target
#
add_custom_target(check_prec COMMAND ${CMAKE_CTEST_COMMAND}
                     DEPENDS ${PrecTests})
