/** \file ParametersFileDefs.cpp
 *  \brief Source of the definitions and names of the parameters file
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/ASCII/ParametersFileDefs.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string ParametersFileDefs::HEADER = "ParametersFile";

   const std::string ParametersFileDefs::VERSION = "1.0";

   const std::string ParametersFileDefs::BASENAME = "parameters";

   const std::string ParametersFileDefs::EXTENSION = ".cfg";

   const std::string ParametersFileDefs::TRUNCXML = "truncation";

   const std::string ParametersFileDefs::TRUNCNXML = "radial";

   const std::string ParametersFileDefs::TRUNCLXML = "latitudinal";

   const std::string ParametersFileDefs::TRUNCMXML = "longitudinal";

   const std::string ParametersFileDefs::TRUNCMPXML = "symmetry";

   const std::string ParametersFileDefs::TRUNCNCOREXML = "cpus";

   const std::string ParametersFileDefs::PHYSICALXML = "physical";

   const std::string ParametersFileDefs::PHYSEKMANXML = "ekman";

   const std::string ParametersFileDefs::PHYSROBERTSXML = "roberts";

   const std::string ParametersFileDefs::PHYSRAYLEIGHXML = "rayleigh";

   const std::string ParametersFileDefs::PHYSROSBYXML = "rosby";

   const std::string ParametersFileDefs::TIMEXML = "timestepping";

   const std::string ParametersFileDefs::TIMETIMEXML = "time";

   const std::string ParametersFileDefs::TIMETSTEPXML = "timestep";

   const std::string ParametersFileDefs::RUNXML = "run";

   const std::string ParametersFileDefs::RUNMAXTSTEPXML = "maxtstep";

   const std::string ParametersFileDefs::RUNARATEXML = "arate";

   const std::string ParametersFileDefs::RUNSRATEXML = "srate";

   const std::string ParametersFileDefs::RUNWALLXML = "walltime";

   ParametersFileDefs::ParametersFileDefs()
   {
   }

}
