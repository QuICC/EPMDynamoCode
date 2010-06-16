/** \file CSCSFileDefs.hpp
 *  \brief Definitions and names for the CSCS visualisation format
 */

#ifndef CSCSFILEDEFS_HPP
#define CSCSFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names for the CSCS visualisation format
    */
   class CSCSFileDefs
   {
      public:
         /**
          * @brief HEADER part for CSCS file
          */
         static const std::string   GRIDTAG;

         /**
          * @brief HEADER part for CSCS file
          */
         static const std::string   RADIALAXISTAG;

         /**
          * @brief HEADER part for CSCS file
          */
         static const std::string   THETAAXISTAG;

         /**
          * @brief HEADER part for CSCS file
          */
         static const std::string   PHIAXISTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   RADIALTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   THETATAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   PHITAG;

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
          * @brief HEADER part for CSCS file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for CSCS file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of CSCS file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of CSCS file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~CSCSFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         CSCSFileDefs();
   };
}

#endif // CSCSFILEDEFS_HPP
