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
          * @brief EXTENSION of parameters file
          */
         static const std::string   TRUNCXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TRUNCNXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TRUNCLXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TRUNCMXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TRUNCMPXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TRUNCNCOREXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   PHYSICALXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   PHYSEKMANXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   PHYSROBERTSXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   PHYSRAYLEIGHXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   PHYSROSSBYXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TIMEXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TIMETIMEXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   TIMETSTEPXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   RUNXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   RUNMAXTSTEPXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   RUNARATEXML;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   RUNSRATEXML;

         /**
          * @brief EXTENSION of parameters file
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
