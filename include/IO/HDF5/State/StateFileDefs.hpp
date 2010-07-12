/** \file StateFileDefs.hpp
 *  \brief Definitions and names use by the state file readers/writers
 */

#ifndef STATEFILEDEFS_HPP
#define STATEFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names use by the state file readers/writers
    */
   class StateFileDefs
   {
      public:

         /**
          * @brief HEADER part for State file
          */
         static const std::string   PHYSICALTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   EKMANTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   ROBERTSTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   RAYLEIGHTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   ROSSBYTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   RUNTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   RUNTIMETAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   RUNSTEPTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   TOROIDALTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   POLOIDALTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   CODENSITYTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   MAGNETICTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   VELOCITYTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for State file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of State file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of State file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~StateFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         StateFileDefs();
   };
}

#endif // STATEFILEDEFS_HPP
