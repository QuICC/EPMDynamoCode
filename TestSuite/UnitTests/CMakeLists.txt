#
# Add Associated Legendre unit test
#
set(UnitTests ${UnitTests} unit_AssociatedLegendre)

#
# Add CSCS file unit test
#
set(UnitTests ${UnitTests} unit_CSCSFile)

#
# Add imposed field unit test
#
set(UnitTests ${UnitTests} unit_ImposedField)

#
# Add parallelisers unit test
#
set(UnitTests ${UnitTests} unit_Parallelisers)

#
# Add ParametersFile unit test
#
set(UnitTests ${UnitTests} unit_ParametersFile)

#
# Add Source file unit test
#
set(UnitTests ${UnitTests} unit_SourceFile)

#
# Add state file unit test
#
set(UnitTests ${UnitTests} unit_StateFile)

#
# Add truncation unit test
#
set(UnitTests ${UnitTests} unit_Truncation)

#
# Add Worland unit test
#
set(UnitTests ${UnitTests} unit_Worland)

#
# Create all the tests
#
foreach(UnitTest ${UnitTests})
   set(SrcsList ${UnitTest} ${All_Srcs})
   add_executable(${UnitTest} ${SrcsList})
   add_test(${UnitTest} ${CMAKE_CURRENT_BINARY_DIR}/${UnitTest})
endforeach(UnitTest)

#
# Create test target
#
add_custom_target(check_unit COMMAND ${CMAKE_CTEST_COMMAND}
                     DEPENDS ${UnitTests})
