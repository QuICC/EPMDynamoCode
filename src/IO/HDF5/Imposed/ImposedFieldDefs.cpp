/** \file ImposedFieldDefs.cpp
 *  \brief Source of the definitions and names used by the imposed field readers/writers
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/Imposed/ImposedFieldDefs.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string ImposedFieldDefs::TOROIDALTAG = "Tor";

   const std::string ImposedFieldDefs::POLOIDALTAG = "Pol";

   const std::string ImposedFieldDefs::CODENSITYTAG = "Codensity";

   const std::string ImposedFieldDefs::MAGNETICTAG = "Magnetic";

   const std::string ImposedFieldDefs::VELOCITYTAG = "Velocity";

   const std::string ImposedFieldDefs::HEADER = "ImposedField";

   const std::string ImposedFieldDefs::VERSION = "1.0";

   const std::string ImposedFieldDefs::BASENAME = "imposed_";

   const std::string ImposedFieldDefs::EXTENSION = ".hdf5";

   ImposedFieldDefs::ImposedFieldDefs()
   {
   }
}
