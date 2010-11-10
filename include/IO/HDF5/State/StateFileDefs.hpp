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

         /**
          * @brief Filter id for not using field
          */
         static const int  DontUse = 0;

         /**
          * @brief Filter id for using full field
          */
         static const int  FullField = 1;

         /**
          * @brief Filter id for using only toroidal component of field
          */
         static const int  ToroidalOnly = 2;

         /**
          * @brief Filter id for using only poloidal component of field
          */
         static const int  PoloidalOnly = 3;

      private:
         /**
         * @brief Empty destructor
         */
         StateFileDefs();
   };
}

#endif // STATEFILEDEFS_HPP
