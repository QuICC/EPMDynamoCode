#
# Add Tubular Splitting para test
#
set(ParaTests ${ParaTests} para_TubularSplitting)

#
# Create all the tests
#
foreach(ParaTest ${ParaTests})
   set(SrcsList ${ParaTest} ${All_Srcs})
   add_executable(${ParaTest} ${SrcsList})
   add_test(${ParaTest} ${CMAKE_CURRENT_BINARY_DIR}/${ParaTest})
endforeach(ParaTest)

#
# Create test target
#
add_custom_target(check_para COMMAND ${CMAKE_CTEST_COMMAND}
                     DEPENDS ${ParaTests})
