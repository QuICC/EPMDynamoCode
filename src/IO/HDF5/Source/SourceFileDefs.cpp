/** \file SourceFileDefs.cpp
 *  \brief Source of the definitions and names used by the source readers/writers
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/Source/SourceFileDefs.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string SourceFileDefs::TOROIDALTAG = "Tor";

   const std::string SourceFileDefs::POLOIDALTAG = "Pol";

   const std::string SourceFileDefs::CODENSITYTAG = "Codensity";

   const std::string SourceFileDefs::MAGNETICTAG = "Magnetic";

   const std::string SourceFileDefs::VELOCITYTAG = "Velocity";

   const std::string SourceFileDefs::HEADER = "Source";

   const std::string SourceFileDefs::VERSION = "1.0";

   const std::string SourceFileDefs::BASENAME = "source_";

   const std::string SourceFileDefs::EXTENSION = ".hdf5";

   SourceFileDefs::SourceFileDefs()
   {
   }
}
