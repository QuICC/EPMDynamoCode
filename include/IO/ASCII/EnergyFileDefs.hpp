/** \file EnergyFileDefs.hpp
 *  \brief Definitions and names used by the Energy output files
 */

#ifndef ENERGYFILEDEFS_HPP
#define ENERGYFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the Energy output files
    */
   class EnergyFileDefs
   {
      public:

         /**
          * @brief HEADER part for energy file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for energy file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of energy file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of energy file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~EnergyFileDefs(){};

      private:
         /**
         * @brief Empty destructor
         */
         EnergyFileDefs();
   };
}

#endif // ENERGYFILEDEFS_HPP
