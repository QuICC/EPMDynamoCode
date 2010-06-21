#
# Create DynamoRun 
#
set(RunName DynamoRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create MagnetoConvectionRun 
#
set(RunName MagnetoConvection)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create ThermalConvectionRun 
#
set(RunName ThermalConvectionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})
