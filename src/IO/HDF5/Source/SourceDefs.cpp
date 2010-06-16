/** \file SourceDefs.cpp
 *  \brief Source of the definitions and names used by the source readers/writers
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/SourceDefs.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string SourceDefs::TOROIDALTAG = "Tor";

   const std::string SourceDefs::POLOIDALTAG = "Pol";

   const std::string SourceDefs::CODENSITYTAG = "Codensity";

   const std::string SourceDefs::MAGNETICTAG = "Magnetic";

   const std::string SourceDefs::VELOCITYTAG = "Velocity";

   const std::string SourceDefs::HEADER = "Source";

   const std::string SourceDefs::VERSION = "1.0";

   const std::string SourceDefs::BASENAME = "source_";

   const std::string SourceDefs::EXTENSION = ".hdf5";

   SourceDefs::SourceDefs()
   {
   }
}
