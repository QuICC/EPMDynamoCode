#
# Project information
#
project(EPMDynamo CXX)

#
# Check for minimal version for cmake
#
cmake_minimum_required(VERSION 2.6.2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#----------------------------- CONFIGURABLE ---------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


###################################################
#--------------------- PLATFORM ------------------#
###################################################

#
# Set platform
# Possible platforms are: Core2 (gentoo), Brutus, Rosa 
#
set(EPMDYNAMO_PLATFORM "${EPMDYNAMO_PLATFORM}" CACHE STRING "Possible platforms: Core2, Brutus, Rosa." FORCE)

# Set Core2 platform
if(EPMDYNAMO_PLATFORM STREQUAL "Core2")
   set(EPMDYNAMO_ON_CORE2 ON)
# Set Brutus platform
elseif(EPMDYNAMO_PLATFORM STREQUAL "Brutus")
   set(EPMDYNAMO_ON_BRUTUS ON)
# Set Rosa platform
elseif(EPMDYNAMO_PLATFORM STREQUAL "Rosa")
   set(EPMDYNAMO_ON_ROSA ON)
else(EPMDYNAMO_PLATFORM STREQUAL "Core2")
   message(SEND_ERROR "Unkown platform")
endif(EPMDYNAMO_PLATFORM STREQUAL "Core2")


###################################################
#--------------------- COMPILER ------------------#
###################################################

#
# Choose compiler
# Possible compilers are: GCC, PGI, Intel, Pathscale, Cray, Sun, Tau, TauInstr
# Not all compilers can be used everywhere.
#
set(EPMDYNAMO_COMPILER "${EPMDYNAMO_COMPILER}" CACHE STRING "Compiler to use: GCC, PGI, Intel, Pathscale, Cray, Sun, Tau, TauInstr." FORCE)

# Set GCC compiler
if(EPMDYNAMO_COMPILER STREQUAL "GCC")
   set(EPMDYNAMO_USE_GCC ON)
# Set PGI compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "PGI")
   set(EPMDYNAMO_USE_PGI ON)
# Set Inte; compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Intel")
   set(EPMDYNAMO_USE_INTEL ON)
# Set Pathscale compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Pathscale")
   set(EPMDYNAMO_USE_PATHSCALE ON)
# Set Cray compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Cray")
   set(EPMDYNAMO_USE_CRAY ON)
# Set Sun compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Sun")
   set(EPMDYNAMO_USE_SUN ON)
# Set Tau compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Tau")
   set(EPMDYNAMO_USE_TAU ON)
# Set Tau compiler with compiler instrumentation
elseif(EPMDYNAMO_COMPILER STREQUAL "TauInstr")
   set(EPMDYNAMO_USE_TAUINSTR ON)
# Unknow compiler
else(EPMDYNAMO_COMPILER STREQUAL "GCC")
   message(SEND_ERROR "Unkown compiler")
endif(EPMDYNAMO_COMPILER STREQUAL "GCC")


###################################################
#---------------- SHARED POINTER -----------------#
###################################################

#
# Choose shared_ptr implementation
# Possible implementations are: Boost, TR1, c++0x
#
set(EPMDYNAMO_SHAREDPTR_IMPLEMENTATION "${EPMDYNAMO_SHAREDPTR_IMPLEMENTATION}" CACHE STRING "Implementations: Boost, TR1, c++0x." FORCE)
MARK_AS_ADVANCED(EPMDYNAMO_SHAREDPTR_IMPLEMENTATION)

# Use Boost shared_ptr implementation
if(EPMDYNAMO_SHAREDPTR_IMPLEMENTATION STREQUAL "Boost")
   add_definitions("-DEPMDYNAMO_BOOST_SMARTPTR")
# Use TR1 shared_ptr implementation
elseif(EPMDYNAMO_SHAREDPTR_IMPLEMENTATION STREQUAL "TR1")
   add_definitions("-DEPMDYNAMO_TR1_SMARTPTR")
# Use C++0x shared_ptr implementation
elseif(EPMDYNAMO_SHAREDPTR_IMPLEMENTATION STREQUAL "c++0x")
   add_definitions("-DEPMDYNAMO_CXX0X_SMARTPTR")
   set(EPMDYNAMO_USE_CXX0X ON)
# Use Boost shared_ptr implementation by default
else(EPMDYNAMO_SHAREDPTR_IMPLEMENTATION STREQUAL "Boost")
   add_definitions("-DEPMDYNAMO_BOOST_SMARTPTR")
endif(EPMDYNAMO_SHAREDPTR_IMPLEMENTATION STREQUAL "Boost")


###################################################
#------------- MPI PARALLELISATION ---------------#
###################################################

#
# Introduce MPI step in Spherical Harmonics transform ?
#
option(EPMDYNAMO_SPLIT_SH "Use MPI spherical harmonics transform version?" OFF)
if(EPMDYNAMO_SPLIT_SH)
   add_definitions("-DEPMDYNAMO_SPLIT_SH")
endif(EPMDYNAMO_SPLIT_SH)

#
# Introduce MPI step in Radial transform ?
#
option(EPMDYNAMO_SPLIT_RADIAL "Use MPI radial transform version?" OFF)
if(EPMDYNAMO_SPLIT_RADIAL)
   add_definitions("-DEPMDYNAMO_SPLIT_RADIAL")
endif(EPMDYNAMO_SPLIT_RADIAL)

#
# Spherical harmonics splitting is only possible if radial transform is also
# split (for the moment)
#
#if(EPMDYNAMO_SPLIT_SH AND NOT EPMDYNAMO_SPLIT_RADIAL)
#   message(SEND_ERROR "Only parallelising the SH transform is not supported yet")
#endif(EPMDYNAMO_SPLIT_SH AND NOT EPMDYNAMO_SPLIT_RADIAL)

#
# If any of the previous settings is on activate MPI flag
#
if(EPMDYNAMO_SPLIT_SH OR EPMDYNAMO_SPLIT_RADIAL)
   set(EPMDYNAMO_MPI ON)
   add_definitions("-DEPMDYNAMO_MPI")
endif(EPMDYNAMO_SPLIT_SH OR EPMDYNAMO_SPLIT_RADIAL)

#
# Use grouped communication in Spherical Harmonics transform ?
#
option(EPMDYNAMO_SH_GROUPEDCOMM "Use grouped MPI communication for spherical harmonics transform" OFF)
if(EPMDYNAMO_SH_GROUPEDCOMM)
   if(EPMDYNAMO_SPLIT_SH)
      add_definitions("-DEPMDYNAMO_SH_GROUPEDCOMM")
   endif(EPMDYNAMO_SPLIT_SH)
endif(EPMDYNAMO_SH_GROUPEDCOMM)

#
# Use grouped communication in Radial transform ?
#
option(EPMDYNAMO_RADIAL_GROUPEDCOMM "Use grouped MPI communication for radial transform?" OFF)
if(EPMDYNAMO_RADIAL_GROUPEDCOMM)
   if(EPMDYNAMO_SPLIT_RADIAL)
      add_definitions("-DEPMDYNAMO_RADIAL_GROUPEDCOMM")
   endif(EPMDYNAMO_SPLIT_RADIAL)
endif(EPMDYNAMO_RADIAL_GROUPEDCOMM)


###################################################
#----------------- MEMORY USAGE ------------------#
###################################################

#
# Choose the type of memory usage the code is allowed to do.
# Possible options are: BigMem, LowMem
#
set(EPMDYNAMO_MEMORY_USAGE "${EPMDYNAMO_MEMORY_USAGE}" CACHE STRING "Possible types of memory usage: BigMem, LowMem." FORCE)

# Set Core2 platform
if(EPMDYNAMO_MEMORY_USAGE STREQUAL "BigMem")
   add_definitions("-DEPMDYNAMO_BIGMEM")
elseif(EPMDYNAMO_MEMORY_USAGE STREQUAL "LowMem")
   add_definitions("-DEPMDYNAMO_LOWMEM")
else(EPMDYNAMO_MEMORY_USAGE STREQUAL "BigMem")
   add_definitions("-DEPMDYNAMO_BIGMEM")
endif(EPMDYNAMO_MEMORY_USAGE STREQUAL "BigMem")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#--------------------------- AUTOMATIC SETUP --------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


###################################################
#--------------------- COMPILER ------------------#
###################################################

#
# General compiler Settings
#
set(CMAKE_CXX_FLAGS_RELEASE "-DEIGEN_NO_DEBUG -DNDEBUG" CACHE STRING "" FORCE)

#
# General GCC compiler settings
#
if(EPMDYNAMO_USE_GCC)
   if(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-O2 -std=c++0x" CACHE STRING "" FORCE)
   else(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-O2" CACHE STRING "" FORCE)
   endif(EPMDYNAMO_USE_CXX0X)
   set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline -funroll-loops -finline-functions" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_GCC)

#
# General PGI compiler Settings
#
if(EPMDYNAMO_USE_PGI)
   #set(CMAKE_CXX_FLAGS "-fastsse")
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Minline" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_PGI)

#
# General Intel compiler settings
#
if(EPMDYNAMO_USE_INTEL)
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_INTEL)

#
# General Pathscale compiler settings
#
if(EPMDYNAMO_USE_PATHSCALE)
   set(CMAKE_CXX_FLAGS "-OPT:Ofast" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_PATHSCALE)

#
# General Cray compiler settings
#
if(EPMDYNAMO_USE_CRAY)
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_CRAY)

#
# General Sun compiler settings
#
if(EPMDYNAMO_USE_SUN)
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_SUN)

#
# General Tau compiler settings
#
if(EPMDYNAMO_USE_TAU)
   if(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-O2 -std=c++0x" CACHE STRING "" FORCE)
   else(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-O2" CACHE STRING "" FORCE)
   endif(EPMDYNAMO_USE_CXX0X)
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline -funroll-loops -finline-functions" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_TAU)

#
# General Tau compiler settings
#
if(EPMDYNAMO_USE_TAUINSTR)
   if(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-tau_options=-optCompInst -O2 -std=c++0x" CACHE STRING "" FORCE)
   else(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-tau_options=-optCompInst -O2" CACHE STRING "" FORCE)
   endif(EPMDYNAMO_USE_CXX0X)
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline -funroll-loops -finline-functions" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_TAUINSTR)

#
# Add libraries to link with
#
if(EPMDYNAMO_ON_CORE2)
   # give error if unavailable compiler is chosen
   if(EPMDYNAMO_USE_PGI)
      message(SEND_ERROR "PGI compiler not available")
   elseif(EPMDYNAMO_USE_INTEL)
      message(SEND_ERROR "Intel compiler not available")
   elseif(EPMDYNAMO_USE_PATHSCALE)
      message(SEND_ERROR "Pathscale compiler not available")
   elseif(EPMDYNAMO_USE_CRAY)
      message(SEND_ERROR "Cray compiler not available")
   elseif(EPMDYNAMO_USE_SUN)
      message(SEND_ERROR "Sun compiler not available")
   elseif(EPMDYNAMO_USE_TAU)
      message(SEND_ERROR "Tau compiler not available")
   elseif(EPMDYNAMO_USE_TAUINSTR)
      message(SEND_ERROR "Tau compiler not available")
   endif(EPMDYNAMO_USE_PGI)

   # Set MPI Wrapper if mpi  is set
   if(EPMDYNAMO_MPI)
      set(CMAKE_CXX_COMPILER "mpicxx")
   else(EPMDYNAMO_MPI)
      if(EPMDYNAMO_USE_GCC)
         set(CMAKE_CXX_COMPILER "g++")
      endif(EPMDYNAMO_USE_GCC)
   endif(EPMDYNAMO_MPI)

   if(EPMDYNAMO_USE_GCC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3 -march=core2" CACHE STRING "" FORCE)
   endif(EPMDYNAMO_USE_GCC)

   link_libraries(lapack)
   link_libraries(fftw3)
   link_libraries(rt)
   link_libraries(hdf5)

endif(EPMDYNAMO_ON_CORE2)

if(EPMDYNAMO_ON_BRUTUS)
   # give error if unavailable compiler is chosen
   if(EPMDYNAMO_USE_PGI)
      message(SEND_ERROR "PGI compiler not tested")
   elseif(EPMDYNAMO_USE_PATHSCALE)
      message(SEND_ERROR "Pathscale compiler not available")
   elseif(EPMDYNAMO_USE_CRAY)
      message(SEND_ERROR "Cray compiler not available")
   elseif(EPMDYNAMO_USE_TAU)
      message(SEND_ERROR "Tau compiler not available")
   elseif(EPMDYNAMO_USE_TAUINSTR)
      message(SEND_ERROR "Tau compiler not available")
   endif(EPMDYNAMO_USE_PGI)

   if(EPMDYNAMO_MPI)
      set(CMAKE_CXX_COMPILER "mpic++")
      if(EPMDYNAMO_USE_SUN)
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
      endif(EPMDYNAMO_USE_SUN)
   else(EPMDYNAMO_MPI)
      if(EPMDYNAMO_USE_GCC)
         set(CMAKE_CXX_COMPILER "g++")
      elseif(EPMDYNAMO_USE_INTEL)
         set(CMAKE_CXX_COMPILER "ifcp")
      elseif(EPMDYNAMO_USE_SUN)
         set(CMAKE_CXX_COMPILER "sunCC")
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
      endif(EPMDYNAMO_USE_GCC)
   endif(EPMDYNAMO_MPI)

   link_libraries(lapack)
   link_libraries(fftw3)
   link_libraries(rt)
   link_libraries(hdf5)
   link_libraries(z)

   # Set special link libraries for GCC
   if(EPMDYNAMO_USE_GCC)
	   include_directories("/cluster/home/erdw/phmarti/myProgs/bin/include")
	   link_directories("/cluster/home/erdw/phmarti/myProgs/bin/lib")
   endif(EPMDYNAMO_USE_GCC)

   # Set special link libraries for Intel
   if(EPMDYNAMO_USE_INTEL)
	   include_directories("/cluster/apps/fftw/x86_64/intel/include")
	   link_directories("/cluster/apps/fftw/x86_64/intel/lib")
   endif(EPMDYNAMO_USE_INTEL)
endif(EPMDYNAMO_ON_BRUTUS)

if(EPMDYNAMO_ON_ROSA)
   # give error if unavailable compiler is chosen
   if(EPMDYNAMO_USE_INTEL)
      message(SEND_ERROR "Intel compiler not  available")
   elseif(EPMDYNAMO_USE_CRAY)
      message(SEND_ERROR "Cray compiler not working")
   elseif(EPMDYNAMO_USE_SUN)
      message(SEND_ERROR "Sun compiler not  available")
   endif(EPMDYNAMO_USE_INTEL)

   if(EPMDYNAMO_MPI)
      set(CMAKE_CXX_COMPILER "CC")
      if(EPMDYNAMO_USE_PGI)
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
      endif(EPMDYNAMO_USE_PGI)
   else(EPMDYNAMO_MPI)
      if(EPMDYNAMO_USE_GCC)
         set(CMAKE_CXX_COMPILER "CC")
      elseif(EPMDYNAMO_USE_PATHSCALE)
         set(CMAKE_CXX_COMPILER "CC")
      elseif(EPMDYNAMO_USE_PGI)
         set(CMAKE_CXX_COMPILER "CC")
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
      endif(EPMDYNAMO_USE_GCC)
   endif(EPMDYNAMO_MPI)

   if(EPMDYNAMO_USE_TAU)
      set(CMAKE_CXX_COMPILER "tau_cxx.sh")
      # Set Boost definitions
	   add_definitions($ENV{BOOST_POST_COMPILE_OPTS})
	   add_definitions("-DMPICH_IGNORE_CXX_SEEK")
   endif(EPMDYNAMO_USE_TAU)

   if(EPMDYNAMO_USE_TAUINSTR)
      set(CMAKE_CXX_COMPILER "tau_cxx.sh")
      # Set Boost definitions
	   add_definitions($ENV{BOOST_POST_COMPILE_OPTS})
	   add_definitions("-DMPICH_IGNORE_CXX_SEEK")
   endif(EPMDYNAMO_USE_TAUINSTR)

   # Set HDF5 definitions
	add_definitions($ENV{HDF5_INCLUDE_OPTS})
	link_directories($ENV{HDF5_DIR}/lib)

   link_libraries(hdf5)
   link_libraries(z)

endif(EPMDYNAMO_ON_ROSA)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#----------------------- CMAKE TESTING STRUCTURE ----------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#
# Enable testing
#
option(EPMDYNAMO_TESTING "Create the CTest tests?" OFF)
if(EPMDYNAMO_TESTING)
   enable_testing()
endif(EPMDYNAMO_TESTING)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#------------------------- DIRECTORY STRUCTURE ------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#
# General variables
# 
# Some directories
set(IncludeDir ${CMAKE_SOURCE_DIR}/include)
set(SrcDir ${CMAKE_SOURCE_DIR}/src)
set(ExternalDir ${CMAKE_SOURCE_DIR}/External)
set(TestSuiteDir ${CMAKE_SOURCE_DIR}/TestSuite)
set(ProgramDir ${CMAKE_SOURCE_DIR}/Programs)
set(EigenDir ${ExternalDir}/eigen2)
set(XMLDir ${ExternalDir}/rapidxml)
set(IFaceDir ${ExternalDir}/Interfaces)
set(CodeDocDir ${CMAKE_BINARY_DIR}/CodeDoc)

#
# Set General include directories
#
include_directories(${IncludeDir})
include_directories(${EigenDir})
include_directories(${XMLDir})
include_directories(${IFaceDir})

#
# Create special rule for doxygen documentation
#
configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.cmake.in
                          ${CMAKE_BINARY_DIR}/Doxyfile)

set(DOXYGEN_INPUT ${CMAKE_BINARY_DIR}/Doxyfile)
set(DOXYGEN_OUTPUT ${CodeDocDir}/html/index.html)
set(DOXYGEN_EXECUTABLE doxygen)

add_custom_command(
   OUTPUT ${DOXYGEN_OUTPUT}
   COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
   COMMAND ${CMAKE_COMMAND} -E echo "Done."
   DEPENDS ${DOXYGEN_INPUT}
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(CodeDoc
   COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
   COMMAND ${CMAKE_COMMAND} -E echo "Done."
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
   )

add_custom_target(cleanCodeDoc
   COMMAND ${CMAKE_COMMAND} -E echo_append "Erasing API Documentation..."
   COMMAND rm -rf ${CodeDocDir}/latex
   COMMAND rm -rf ${CodeDocDir}/html
   COMMAND ${CMAKE_COMMAND} -E echo "Done."
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
   )


#
# Generate Source file list
#
set(SrcSubDirs
      BoundaryConditions
      DataManipulators
      Domain
      Equations
      General
      GeneralFields
      GeneralScalars
      GridGenerators
      IO
      Operators
      Parallelisers
      Polynomials
      Simulations
      Timestepping
      Transforms)

foreach(SrcSubDir ${SrcSubDirs})
   set(SrcSubDirPath ${SrcDir}/${SrcSubDir})
   include(${SrcSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_Srcs})
      list(APPEND All_Srcs ${SrcSubDirPath}/${source})
   endforeach(source)
endforeach(SrcSubDir)

#
# Add project subdirectories
#
add_subdirectory(${SrcDir})
add_subdirectory(${ProgramDir} EXCLUDE_FROM_ALL)
add_subdirectory(${TestSuiteDir} EXCLUDE_FROM_ALL)

