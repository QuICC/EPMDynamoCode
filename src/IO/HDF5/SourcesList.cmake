#
# Create sources list
#
set(${SrcSubDir}_${SrcSubSubDir}_Srcs
      HDF5File.cpp
      HDF5FileDefs.cpp
      HDF5NWriter.cpp
      HDF5Reader.cpp
      HDF5TypesConv.cpp
      HDF5Writer.cpp
      SpectralFileDefs.cpp
      SpectralHDF5NWriter.cpp
      SpectralHDF5Reader.cpp
      )

#
# include source files from subsubdirectories
#
set(SrcSubSubSubDirs
      CSCS
      Imposed
      Source
      State
      )

foreach(SrcSubSubSubDir ${SrcSubSubSubDirs})
   set(SrcSubSubSubDirPath ${SrcSubSubDirPath}/${SrcSubSubSubDir})
   include(${SrcSubSubSubDirPath}/SourcesList.cmake)
   foreach(source ${${SrcSubDir}_${SrcSubSubDir}_${SrcSubSubSubDir}_Srcs})
      list(APPEND ${SrcSubDir}_${SrcSubSubDir}_Srcs ${SrcSubSubSubDir}/${source})
   endforeach(source)
endforeach(SrcSubSubSubDir)
