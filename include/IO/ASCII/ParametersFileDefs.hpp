/** \file ParametersFileDefs.hpp
 *  \brief Definitions and names used by the parameters file
 */

#ifndef PARAMETERSFILEDEFS_HPP
#define PARAMETERSFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the parameters file
    */
   class ParametersFileDefs
   {
      public:

         /**
          * @brief HEADER part for parameters file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for parameters file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of parameters file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   EXTENSION;

         /**
          * @brief Truncation part of parameters file
          */
         static const std::string   TRUNCXML;

         /**
          * @brief Radial truncation node of parameters file
          */
         static const std::string   TRUNCNXML;

         /**
          * @brief Harmonic degree truncation node of parameters file
          */
         static const std::string   TRUNCLXML;

         /**
          * @brief Harmonic order truncation node of parameters file
          */
         static const std::string   TRUNCMXML;

         /**
          * @brief Assumed symmetry node of parameters file
          */
         static const std::string   TRUNCMPXML;

         /**
          * @brief Number of cpu node of parameters file
          */
         static const std::string   TRUNCNCOREXML;

         /**
          * @brief Physical parameters node of parameters file
          */
         static const std::string   PHYSICALXML;

         /**
          * @brief Ekman number node of parameters file
          */
         static const std::string   PHYSEKMANXML;

         /**
          * @brief Elsasser number node of parameters file
          */
         static const std::string   PHYSELSASSERXML;

         /**
          * @brief Thermal Ekman number node of parameters file
          */
         static const std::string   PHYSTHERMALEKMANXML;

         /**
          * @brief Magnetic Ekman number node of parameters file
          */
         static const std::string   PHYSMAGNETICEKMANXML;

         /**
          * @brief Magnetic Prandtl number node of parameters file
          */
         static const std::string   PHYSMAGNETICPRANDTLXML;

         /**
          * @brief Prandtl number node of parameters file
          */
         static const std::string   PHYSPRANDTLXML;

         /**
          * @brief Roberts number node of parameters file
          */
         static const std::string   PHYSROBERTSXML;

         /**
          * @brief Rayleigh number node of parameters file
          */
         static const std::string   PHYSRAYLEIGHXML;

         /**
          * @brief Rossby number node of parameters file
          */
         static const std::string   PHYSROSSBYXML;

         /**
          * @brief Boundary condition node of parameters file
          */
         static const std::string   BOUNDARYXML;

         /**
          * @brief Codensity boundary condition of parameters file
          */
         static const std::string   BCCODXML;

         /**
          * @brief Velocity boundary condition of parameters file
          */
         static const std::string   BCVELXML;

         /**
          * @brief Magnetic boundary condition of parameters file
          */
         static const std::string   BCMAGXML;

         /**
          * @brief Timestepping node of parameters file
          */
         static const std::string   TIMEXML;

         /**
          * @brief Starting time node of parameters file
          */
         static const std::string   TIMETIMEXML;

         /**
          * @brief Timestep node of parameters file
          */
         static const std::string   TIMETSTEPXML;

         /**
          * @brief Runtime node of parameters file
          */
         static const std::string   RUNXML;

         /**
          * @brief Maximum integration time node of parameters file
          */
         static const std::string   RUNMAXTIMEXML;

         /**
          * @brief Maximum number of timesteps node of parameters file
          */
         static const std::string   RUNMAXTSTEPXML;

         /**
          * @brief ASCII diagnostics save rate node of parameters file
          */
         static const std::string   RUNARATEXML;

         /**
          * @brief State file save rate node of parameters file
          */
         static const std::string   RUNSRATEXML;

         /**
          * @brief Wall time node of parameters file
          */
         static const std::string   RUNWALLXML;

         /**
         * @brief Destructor
         */
         virtual ~ParametersFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         ParametersFileDefs();
   };
}

#endif // PARAMETERSFILEDEFS_HPP
