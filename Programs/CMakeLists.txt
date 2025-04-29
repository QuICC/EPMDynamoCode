#
# Create all the simulation targets
#
file(GLOB SimImplementations
   ${IncludeDir}/Simulations/Implementations/*Simulation.hpp)
foreach(SimFile ${SimImplementations})
   get_filename_component(SimName ${SimFile} NAME_WE)
   set(SrcsList SimulationRun ${SrcDir}/Simulations/Implementations/${SimName}.cpp ${All_Srcs})
   add_executable(${SimName} ${SrcsList})
   set_target_properties(${SimName} PROPERTIES OUTPUT_NAME ${SimName}
      COMPILE_FLAGS "-D IMPLEMENTATION=${SimName}")
endforeach(SimFile)

#
# Create GenerateCSCS 
#
set(RunName GenerateCSCS)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create GenerateImposed 
#
set(RunName GenerateImposed)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create GenerateSource 
#
set(RunName GenerateSource)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create GenerateState 
#
set(RunName GenerateState)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Compute field integral 
#
set(RunName ComputeIntegral)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# analyse stored state 
#
set(RunName AnalyseState)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})
