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
# Possible platforms are: Core2 (Gentoo), Brutus, Rosa 
#
set(EPMDYNAMO_PLATFORM "${EPMDYNAMO_PLATFORM}" CACHE STRING "Possible platforms: Core2, Brutus, Rosa." FORCE)

# Set Core2 platform
if(EPMDYNAMO_PLATFORM STREQUAL "Core2")
   set(EPMDYNAMO_ON_CORE2 ON)
   message(STATUS "***** Setup system for \"Core2\" platform *****")
# Set Brutus platform
elseif(EPMDYNAMO_PLATFORM STREQUAL "Brutus")
   set(EPMDYNAMO_ON_BRUTUS ON)
   message(STATUS "***** Setup system for \"Brutus\" platform *****")
# Set Rosa platform
elseif(EPMDYNAMO_PLATFORM STREQUAL "Rosa")
   set(EPMDYNAMO_ON_ROSA ON)
   message(STATUS "***** Setup system for \"Rosa\" platform *****")
# Error message if unknown platform is set
else(EPMDYNAMO_PLATFORM STREQUAL "Core2")
   message(SEND_ERROR "Unkown platform")
endif(EPMDYNAMO_PLATFORM STREQUAL "Core2")


###################################################
#--------------------- COMPILER ------------------#
###################################################

#
# Choose compiler
# Possible compilers are: GCC, PGI, Intel, Pathscale, Cray, Sun, TauGCC, TauGCC_intrument, ScalGCC
# Not all compilers can be used everywhere.
#
set(EPMDYNAMO_COMPILER "${EPMDYNAMO_COMPILER}" CACHE STRING "Compiler to use: GCC, PGI, Intel, Pathscale, Cray, Sun, TauGCC, TauGCC_instrument, ScalGCC." FORCE)

# Set GCC compiler
if(EPMDYNAMO_COMPILER STREQUAL "GCC")
   set(EPMDYNAMO_USE_GCC ON)
   message(STATUS "***** Setup system with \"GCC\" compiler *****")
# Set PGI compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "PGI")
   set(EPMDYNAMO_USE_PGI ON)
   message(STATUS "***** Setup system with \"PGI\" compiler *****")
# Set Inte; compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Intel")
   set(EPMDYNAMO_USE_INTEL ON)
   message(STATUS "***** Setup system with \"Intel\" compiler *****")
# Set Pathscale compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Pathscale")
   set(EPMDYNAMO_USE_PATHSCALE ON)
   message(STATUS "***** Setup system with \"Pathscale\" compiler *****")
# Set Cray compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Cray")
   set(EPMDYNAMO_USE_CRAY ON)
   message(STATUS "***** Setup system with \"Cray\" compiler *****")
# Set Sun compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "Sun")
   set(EPMDYNAMO_USE_SUN ON)
   message(STATUS "***** Setup system with \"Sun\" compiler *****")
# Set Tau compiler based on GCC
elseif(EPMDYNAMO_COMPILER STREQUAL "TauGCC")
   set(EPMDYNAMO_USE_TAU ON)
   set(EPMDYNAMO_USE_GCC ON)
   message(STATUS "***** Setup system with \"Tau with GCC\" compiler *****")
# Set Tau compiler based on GCC with compiler instrumentation setup
elseif(EPMDYNAMO_COMPILER STREQUAL "TauGCC_instrument")
   set(EPMDYNAMO_USE_TAU ON)
   set(EPMDYNAMO_USE_TAU_INSTRUMENTED ON)
   set(EPMDYNAMO_USE_GCC ON)
   message(STATUS "***** Setup system with \"Tau with GCC instrumentation\" compiler *****")
# Set Scalasca compiler with gnu base compiler
elseif(EPMDYNAMO_COMPILER STREQUAL "ScalGCC")
   set(EPMDYNAMO_USE_GCC ON)
   set(EPMDYNAMO_USE_SCALASCA ON)
   message(STATUS "***** Setup system with \"Scalasca with GCC\" compiler *****")
# Error message for unknow compiler
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
#------------ FFTW PLAN COMPUTATION --------------#
###################################################

#
# Choose the fftw plan computation flag
# Possible flags are: fast, medium, slow
#
set(EPMDYNAMO_FFTW_PLAN_FLAG "${EPMDYNAMO_FFTW_PLAN_FLAG}" CACHE STRING "FFTW flag: fast, medium, slow." FORCE)
MARK_AS_ADVANCED(EPMDYNAMO_FFTW_PLAN_FLAG)

# Use ESTIMATE FFTW plan flag
if(EPMDYNAMO_FFTW_PLAN_FLAG STREQUAL "fast")
   message(STATUS "Setup fast FFTW plan computation (FFTW_ESTIMATE)")
   add_definitions("-DEPMDYNAMO_FFTW_ESTIMATE")
# Use MEASURE FFTW plan flag
elseif(EPMDYNAMO_FFTW_PLAN_FLAG STREQUAL "medium")
   message(STATUS "Setup normal FFTW plan computation (FFTW_MEASURE)")
   add_definitions("-DEPMDYNAMO_FFTW_MEASURE")
# Use PATIENT FFTW plan flag
elseif(EPMDYNAMO_FFTW_PLAN_FLAG STREQUAL "slow")
   message(STATUS "Setup slow FFTW plan computation (FFTW_PATIENT)")
   add_definitions("-DEPMDYNAMO_FFTW_PATIENT")
# Use ESTIMATE FFTW plan flag by default
else(EPMDYNAMO_FFTW_PLAN_FLAG STREQUAL "fast")
   message(STATUS "Setup fast FFTW plan computation (FFTW_ESTIMATE)")
   add_definitions("-DEPMDYNAMO_FFTW_ESTIMATE")
endif(EPMDYNAMO_FFTW_PLAN_FLAG STREQUAL "fast")


###################################################
#------------- MPI PARALLELISATION ---------------#
###################################################

#
# Introduce MPI step in Spherical Harmonics transform ?
#
option(EPMDYNAMO_SPLIT_SH "Use MPI spherical harmonics transform version?" OFF)
if(EPMDYNAMO_SPLIT_SH)
   message(STATUS "Setting up MPI splitted spherical harmonics transform")
   add_definitions("-DEPMDYNAMO_SPLIT_SH")
endif(EPMDYNAMO_SPLIT_SH)

#
# Introduce MPI step in Radial transform ?
#
option(EPMDYNAMO_SPLIT_RADIAL "Use MPI radial transform version?" OFF)
if(EPMDYNAMO_SPLIT_RADIAL)
   message(STATUS "Setting up MPI splitted radial transform")
   add_definitions("-DEPMDYNAMO_SPLIT_RADIAL")
endif(EPMDYNAMO_SPLIT_RADIAL)

#
# Spherical harmonics splitting is only possible if radial transform is also
# split (for the moment)
#
if(EPMDYNAMO_SPLIT_SH AND NOT EPMDYNAMO_SPLIT_RADIAL)
   message(SEND_ERROR "Only parallelising the SH transform is not supported yet")
endif(EPMDYNAMO_SPLIT_SH AND NOT EPMDYNAMO_SPLIT_RADIAL)

#
# If any of the previous settings is on, activate MPI flag
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
      message(STATUS "Setting up grouped spherical harmonics transform communication")
      add_definitions("-DEPMDYNAMO_SH_GROUPEDCOMM")
   endif(EPMDYNAMO_SPLIT_SH)
endif(EPMDYNAMO_SH_GROUPEDCOMM)

#
# Use grouped communication in Radial transform ?
#
option(EPMDYNAMO_RADIAL_GROUPEDCOMM "Use grouped MPI communication for radial transform?" OFF)
if(EPMDYNAMO_RADIAL_GROUPEDCOMM)
   if(EPMDYNAMO_SPLIT_RADIAL)
      message(STATUS "Setting up grouped radial transform communication")
      add_definitions("-DEPMDYNAMO_RADIAL_GROUPEDCOMM")
   endif(EPMDYNAMO_SPLIT_RADIAL)
endif(EPMDYNAMO_RADIAL_GROUPEDCOMM)


###################################################
#----------------- MEMORY USAGE ------------------#
###################################################

#
# Choose the type of memory usage setup the code is using.
# Possible options are: BigMem, LowMem
#
set(EPMDYNAMO_MEMORY_USAGE "${EPMDYNAMO_MEMORY_USAGE}" CACHE STRING "Possible setups of memory usage: BigMem, LowMem." FORCE)

# Set up big memory usage version
if(EPMDYNAMO_MEMORY_USAGE STREQUAL "BigMem")
   message(STATUS "Setting up big memory usage version")
   add_definitions("-DEPMDYNAMO_BIGMEM")
# Set up low memory usage version
elseif(EPMDYNAMO_MEMORY_USAGE STREQUAL "LowMem")
   add_definitions("-DEPMDYNAMO_LOWMEM")
   message(STATUS "Setting up low memory usage version")
# If unset or unknown value set up big memory usage version
else(EPMDYNAMO_MEMORY_USAGE STREQUAL "BigMem")
   message(STATUS "Memory usage unset. Setting up big memory usage version")
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
   # Set up C++x0 GCC options
   if(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-O2 -std=c++0x" CACHE STRING "" FORCE)
   # Set up standard GCC options
   else(EPMDYNAMO_USE_CXX0X)
      set(CMAKE_CXX_FLAGS "-O2" CACHE STRING "" FORCE)
   endif(EPMDYNAMO_USE_CXX0X)
   # Set debug flags
   set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g" CACHE STRING "" FORCE)
   # Set release flags
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline -funroll-loops -finline-functions" CACHE STRING "" FORCE)
   # Set release flags with debug information
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   # Set profiling flags 
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -pg" CACHE STRING "Flags for doing a profiling" FORCE)
   set(CMAKE_EXE_LINKER_FLAGS_PROFILE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -pg" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE
                     CMAKE_EXE_LINKER_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_GCC)

#
# General PGI compiler Settings
#
if(EPMDYNAMO_USE_PGI)
   # Set basic flags
      #set(CMAKE_CXX_FLAGS "-fastsse")
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   # Set debug flags
   set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "" FORCE)
   # Set release flags
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Minline" CACHE STRING "" FORCE)
   # Set release flags with debug information
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   # Set profiling flags
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_PGI)

#
# General Intel compiler settings
#
if(EPMDYNAMO_USE_INTEL)
   # Set basic flags
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   # Set debug flags
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   # Set release flags
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   # Set release flags with debug information
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   # Set profiling flags
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_INTEL)

#
# General Pathscale compiler settings
#
if(EPMDYNAMO_USE_PATHSCALE)
   # Set basic flags
   set(CMAKE_CXX_FLAGS "-OPT:Ofast" CACHE STRING "" FORCE)
   # Set debug flags
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   # Set release flags
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   # Set release flags with debug information
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   # Set profiling flags
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_PATHSCALE)

#
# General Cray compiler settings
#
if(EPMDYNAMO_USE_CRAY)
   # Set basic flags
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   # Set debug flags
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   # Set release flags
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   # Set release flags with debug information
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   # Set profiling flags
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_CRAY)

#
# General Sun compiler settings
#
if(EPMDYNAMO_USE_SUN)
   # Set basic flags
   set(CMAKE_CXX_FLAGS "" CACHE STRING "" FORCE)
   # Set debug flags
   set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
   # Set release flags
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
   # Set release flags with debug information
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
   # Set profiling flags
   set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE STRING "Flags for doing a profiling" FORCE)
   MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILE)
endif(EPMDYNAMO_USE_SUN)

#
# General Tau compiler settings
#
if(EPMDYNAMO_USE_TAU)
   # Set basic TAU flags
   set(EPMDYNAMO_TAU_FLAGS "")
   # Include TAU flags into basic compiler flags
   set(CMAKE_CXX_FLAGS "${EPMDYNAMO_TAU_FLAGS} ${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE)
endif(EPMDYNAMO_USE_TAU)

#
# General Tau compiler settings
#
if(EPMDYNAMO_USE_TAU_INSTRUMENTED)
   # Set basic TAU flags
   set(EPMDYNAMO_TAU_INSTRUMENTED_FLAGS "-tau_options=-optCompInst")
   # Include TAU compiler instrumented flags into basic compiler flags
   set(CMAKE_CXX_FLAGS "${EPMDYNAMO_TAU_INSTRUMENTED_FLAGS} ${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE)
endif(EPMDYNAMO_USE_TAU_INSTRUMENTED)

#
# General Scalasca compiler settings
#
if(EPMDYNAMO_USE_SCALASCA)
endif(EPMDYNAMO_USE_SCALASCA)

#
# Add libraries to link with
#
if(EPMDYNAMO_ON_CORE2)
   # Error message for unavailable PGI compiler
   if(EPMDYNAMO_USE_PGI)
      message(SEND_ERROR "PGI compiler is not available")
   # Error message for unavailable Intel compiler
   elseif(EPMDYNAMO_USE_INTEL)
      message(SEND_ERROR "Intel compiler is not available")
   # Error message for unavailable Pathscale compiler
   elseif(EPMDYNAMO_USE_PATHSCALE)
      message(SEND_ERROR "Pathscale compiler is not available")
   # Error message for unavailable Cray compiler
   elseif(EPMDYNAMO_USE_CRAY)
      message(SEND_ERROR "Cray compiler is not available")
   # Error message for unavailable Sun compiler
   elseif(EPMDYNAMO_USE_SUN)
      message(SEND_ERROR "Sun compiler is not available")
   # Error message for unavailable TAU system
   elseif(EPMDYNAMO_USE_TAU)
      message(SEND_ERROR "Tau compiler is not available")
   # Error message for unavailable TAU system with compiler instrumentation
   elseif(EPMDYNAMO_USE_TAU_INSTRUMENTED)
      message(SEND_ERROR "Tau compiler is not available")
   # Error message for unavailable Scalasca profiling system
   elseif(EPMDYNAMO_USE_SCALASCA)
      message(SEND_ERROR "Scalasca is not available")
   endif(EPMDYNAMO_USE_PGI)

   # Set MPI Wrapper if required
   if(EPMDYNAMO_MPI)
      set(CMAKE_CXX_COMPILER "mpicxx")
   # Setup serial compilers
   else(EPMDYNAMO_MPI)
      if(EPMDYNAMO_USE_GCC)
         set(CMAKE_CXX_COMPILER "g++")
      endif(EPMDYNAMO_USE_GCC)
   endif(EPMDYNAMO_MPI)

   # Set up platform specific compiler flags
   if(EPMDYNAMO_USE_GCC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3 -march=core2" CACHE STRING "" FORCE)
   endif(EPMDYNAMO_USE_GCC)

   # Link to lapack library
   link_libraries(lapack)
   # Link to FFTW3 library
   link_libraries(fftw3)
   # Link to rt library
   link_libraries(rt)
   # Link to HDF5 library
   link_libraries(hdf5)

endif(EPMDYNAMO_ON_CORE2)

if(EPMDYNAMO_ON_BRUTUS)
   # Error message for unavailable PGI compiler
   if(EPMDYNAMO_USE_PGI)
      message(SEND_ERROR "PGI compiler is not available")
   # Error message for unavailable Pathscale compiler
   elseif(EPMDYNAMO_USE_PATHSCALE)
      message(SEND_ERROR "Pathscale compiler is not available")
   # Error message for unavailable Cray compiler
   elseif(EPMDYNAMO_USE_CRAY)
      message(SEND_ERROR "Cray compiler is not available")
   # Error message for unavailable TAU system
   elseif(EPMDYNAMO_USE_TAU)
      message(SEND_ERROR "Tau compiler is not available")
   # Error message for unavailable TAU system with compiler instrumentation
   elseif(EPMDYNAMO_USE_TAU_INSTRUMENTED)
      message(SEND_ERROR "Tau compiler is not available")
   # Error message for unavailable Scalasca profiling system
   elseif(EPMDYNAMO_USE_SCALASCA)
      message(SEND_ERROR "Scalasca is not available")
   endif(EPMDYNAMO_USE_PGI)

   # Setup MPIwrappers if required
   if(EPMDYNAMO_MPI)
      set(CMAKE_CXX_COMPILER "mpic++")
      if(EPMDYNAMO_USE_SUN)
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
      endif(EPMDYNAMO_USE_SUN)
   # Setup serial compilers
   else(EPMDYNAMO_MPI)
      # Set up GCC compiler
      if(EPMDYNAMO_USE_GCC)
         set(CMAKE_CXX_COMPILER "g++")
      # Set up Intel compiler
      elseif(EPMDYNAMO_USE_INTEL)
         set(CMAKE_CXX_COMPILER "icc")
      # Set up Sun compiler
      elseif(EPMDYNAMO_USE_SUN)
         set(CMAKE_CXX_COMPILER "sunCC")
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
      endif(EPMDYNAMO_USE_GCC)
   endif(EPMDYNAMO_MPI)

   # Link to lapack library
   link_libraries(lapack)
   # Link to FFTW3 library
   link_libraries(fftw3)
   # Link to rt library
   link_libraries(rt)
   # Link to HDF5 library
   link_libraries(hdf5)
   # Link to ZLib library
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
   # Error message for unavailable Intel compiler
   if(EPMDYNAMO_USE_INTEL)
      message(SEND_ERROR "Intel compiler is not available")
   # Error message for unavailable Cray compiler
   elseif(EPMDYNAMO_USE_CRAY)
      message(SEND_ERROR "Cray compiler is not available")
   # Error message for unavailable Sun compiler
   elseif(EPMDYNAMO_USE_SUN)
      message(SEND_ERROR "Sun compiler is not available")
   endif(EPMDYNAMO_USE_INTEL)

   # Change general compilation setup for scalasca
   if(EPMDYNAMO_USE_SCALASCA)
      # Change compilation rule to include scalasca instrumentation
      SET(CMAKE_CXX_COMPILE_OBJECT
         "scalasca -instrument <CMAKE_CXX_COMPILER>  <DEFINES> <FLAGS> -o <OBJECT> -c <SOURCE>")
      # Change executable linkage rule to include scalasca instrumentation
      SET(CMAKE_CXX_LINK_EXECUTABLE
         "scalasca -instrument <CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> <LINK_LIBRARIES>")
   endif(EPMDYNAMO_USE_SCALASCA)

   # Set compiler wrapper independently of chose compiler
   set(CMAKE_CXX_COMPILER "CC")

   # Setup MPI specific variables
   if(EPMDYNAMO_MPI)
      # setup GCC specific variables
      if(EPMDYNAMO_USE_GCC)
         set(HDF5_COMPILER_SPECIFIC_DIR "hdf5-parallel-gnu")
      # setup Pathscale specific variables
      elseif(EPMDYNAMO_USE_PATHSCALE)
         set(HDF5_COMPILER_SPECIFIC_DIR "hdf5-parallel-pathscale")
      # setup PGI specific variables
      elseif(EPMDYNAMO_USE_PGI)
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
         set(HDF5_COMPILER_SPECIFIC_DIR "hdf5-parallel-pgi")
      endif(EPMDYNAMO_USE_GCC)
   # Setup MPI specific variables
   else(EPMDYNAMO_MPI)
      # setup GCC specific variables
      if(EPMDYNAMO_USE_GCC)
         set(HDF5_COMPILER_SPECIFIC_DIR "hdf5-gnu")
      # setup Pathscale specific variables
      elseif(EPMDYNAMO_USE_PATHSCALE)
         set(HDF5_COMPILER_SPECIFIC_DIR "hdf5-pathscale")
      # setup PGI specific variables
      elseif(EPMDYNAMO_USE_PGI)
         set(CMAKE_CXX_FLAGS "-DEIGEN_DONT_ALIGN")
         set(HDF5_COMPILER_SPECIFIC_DIR "hdf5-pgi")
      endif(EPMDYNAMO_USE_GCC)
   endif(EPMDYNAMO_MPI)

   # Set general tau specific settings
   if(EPMDYNAMO_USE_TAU)
      # Set up TAU compiler wrapper
      set(CMAKE_CXX_COMPILER "tau_cxx.sh")
      # add missing boost compilation flags
	   add_definitions($ENV{BOOST_POST_COMPILE_OPTS})
      # add missing FFTW compilation flags
	   add_definitions($ENV{FFTW_INCLUDE_OPTS})
      # add missing FFTW linking flags
	   add_definitions($ENV{FFTW_POST_COMPILE_OPTS})
      # Remove MPICH  SEEK TAU specific multiple definition error
	   add_definitions("-DMPICH_IGNORE_CXX_SEEK")
   endif(EPMDYNAMO_USE_TAU)

   # Set HDF5 library definitions (due to C library for C++ code
	add_definitions(-I$ENV{CRAY_HDF5_DIR}/${HDF5_COMPILER_SPECIFIC_DIR}/include)
	link_directories($ENV{CRAY_HDF5_DIR}/${HDF5_COMPILER_SPECIFIC_DIR}/lib)

   # Link to HDF5 library
   link_libraries(hdf5)
   # Link to ZLib library
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
set(TestSuiteIncDir ${TestSuiteDir}/Tools/include)
set(TestSuiteSrcDir ${TestSuiteDir}/Tools/src)
set(TestSuitePrecIncDir ${TestSuiteDir}/PrecisionTests/Tools/include)
set(TestSuitePrecSrcDir ${TestSuiteDir}/PrecisionTests/Tools/src)
set(ProgramDir ${CMAKE_SOURCE_DIR}/Programs)
set(EigenDir ${ExternalDir}/eigen3)
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
include_directories(${TestSuiteIncDir})
include_directories(${TestSuitePrecIncDir})

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
# Generate general tools Source file list
#
set(ToolsSrcDirPath ${TestSuiteSrcDir})
include(${ToolsSrcDirPath}/SourcesList.cmake)
foreach(source ${${ToolsSrcDir}_Srcs})
   list(APPEND Tools_Srcs ${ToolsSrcDirPath}/${source})
endforeach(source)

#
# Generate PrecisionTests tools Source file list
#
set(ToolsSrcDirPath ${TestSuitePrecSrcDir})
include(${ToolsSrcDirPath}/SourcesList.cmake)
foreach(source ${${ToolsSrcDir}_Srcs})
   list(APPEND Tools_Srcs ${ToolsSrcDirPath}/${source})
endforeach(source)

#
# Add project subdirectories
#
add_subdirectory(${ProgramDir} EXCLUDE_FROM_ALL)
add_subdirectory(${TestSuiteDir} EXCLUDE_FROM_ALL)

