/** \file CSCSFileDefs.cpp
 *  \brief Source of definitions and names used by the CSCS visualisation format
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/Visualisation/CSCSFileDefs.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string CSCSFileDefs::GRIDTAG = "grid";

   const std::string CSCSFileDefs::RADIALAXISTAG = "r axis";

   const std::string CSCSFileDefs::THETAAXISTAG = "t axis";

   const std::string CSCSFileDefs::PHIAXISTAG = "p axis";

   const std::string CSCSFileDefs::RADIALTAG = "r";

   const std::string CSCSFileDefs::THETATAG = "t";

   const std::string CSCSFileDefs::PHITAG = "p";

   const std::string CSCSFileDefs::CODENSITYTAG = "Codensity";

   const std::string CSCSFileDefs::MAGNETICTAG = "Magnetic";

   const std::string CSCSFileDefs::VELOCITYTAG = "Velocity";

   const std::string CSCSFileDefs::HEADER = "CSCSFile";

   const std::string CSCSFileDefs::VERSION = "1.0";

   const std::string CSCSFileDefs::BASENAME = "stateCSCS";

   const std::string CSCSFileDefs::EXTENSION = ".hdf5";

   CSCSFileDefs::CSCSFileDefs()
   {
   }

}
