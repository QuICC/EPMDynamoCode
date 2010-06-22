#
# Create DynamoRun 
#
set(RunName DynamoRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create MagnetoConvectionRun 
#
set(RunName MagnetoConvectionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create ThermalConvectionRun 
#
set(RunName ThermalConvectionRun)
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
