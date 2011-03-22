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
# Create ThermalRotConvectionRun 
#
set(RunName ThermalRotConvectionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create ThermalNLConvectionRun 
#
set(RunName ThermalNLConvectionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create ThermalNLRotConvectionRun 
#
set(RunName ThermalNLRotConvectionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create VelocityDiffusionRun 
#
set(RunName VelocityDiffusionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create VelocityRotDiffusionRun 
#
set(RunName VelocityRotDiffusionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create VelocityNLDiffusionRun 
#
set(RunName VelocityNLDiffusionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create VelocityNLRotDiffusionRun 
#
set(RunName VelocityNLRotDiffusionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create MagneticDiffusionRun 
#
set(RunName MagneticDiffusionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create CodensityDiffusionRun 
#
set(RunName CodensityDiffusionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

#
# Create CodensitySrcDiffusionRun 
#
set(RunName CodensitySrcDiffusionRun)
set(SrcsList ${RunName} ${All_Srcs})
add_executable(${RunName} ${SrcsList})

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
