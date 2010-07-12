/** \file StateFileDefs.cpp
 *  \brief Source of the definitions and names used by the state file readers/writers
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/State/StateFileDefs.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string StateFileDefs::PHYSICALTAG = "PhysicalParameters";

   const std::string StateFileDefs::EKMANTAG = "E";

   const std::string StateFileDefs::ROBERTSTAG = "q";

   const std::string StateFileDefs::RAYLEIGHTAG = "Ra";

   const std::string StateFileDefs::ROSSBYTAG = "Ro";

   const std::string StateFileDefs::RUNTAG = "RunParameters";

   const std::string StateFileDefs::RUNTIMETAG = "Time";

   const std::string StateFileDefs::RUNSTEPTAG = "Step";

   const std::string StateFileDefs::TOROIDALTAG = "Tor";

   const std::string StateFileDefs::POLOIDALTAG = "Pol";

   const std::string StateFileDefs::CODENSITYTAG = "Codensity";

   const std::string StateFileDefs::MAGNETICTAG = "Magnetic";

   const std::string StateFileDefs::VELOCITYTAG = "Velocity";

   const std::string StateFileDefs::HEADER = "StateFile";

   const std::string StateFileDefs::VERSION = "1.0";

   const std::string StateFileDefs::BASENAME = "state";

   const std::string StateFileDefs::EXTENSION = ".hdf5";

   StateFileDefs::StateFileDefs()
   {
   }
}
