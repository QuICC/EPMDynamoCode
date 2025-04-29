/** \file SpectralFileDefs.cpp
 *  \brief Source of the definitions and name used by the spectral data files
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/SpectralFileDefs.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string SpectralFileDefs::TRUNCATIONTAG = "Truncation";

   const std::string SpectralFileDefs::TRUNCNTAG = "N";

   const std::string SpectralFileDefs::TRUNCLTAG = "L";

   const std::string SpectralFileDefs::TRUNCMTAG = "M";

   const std::string SpectralFileDefs::TRUNCMPTAG = "Mp";

   SpectralFileDefs::SpectralFileDefs()
   {
   }
}
