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
          * @brief Grid tag name for CSCS file
          */
         static const std::string   GRIDTAG;

         /**
          * @brief Radial axis tag name for CSCS file
          */
         static const std::string   RADIALAXISTAG;

         /**
          * @brief Theta axis tag name for CSCS file
          */
         static const std::string   THETAAXISTAG;

         /**
          * @brief Phi axis tag name for CSCS file
          */
         static const std::string   PHIAXISTAG;

         /**
          * @brief Radial component tag name for CSCS file
          */
         static const std::string   RADIALTAG;

         /**
          * @brief Theta component tag name for CSCS file
          */
         static const std::string   THETATAG;

         /**
          * @brief Phi component tag name for CSCS file
          */
         static const std::string   PHITAG;

         /**
          * @brief Codensity tag name for CSCS file
          */
         static const std::string   CODENSITYTAG;

         /**
          * @brief Magnetic tag name for CSCS file
          */
         static const std::string   MAGNETICTAG;

         /**
          * @brief Velocity tag name for CSCS file
          */
         static const std::string   VELOCITYTAG;

         /**
          * @brief Vorticity tag name for CSCS file
          */
         static const std::string   VORTICITYTAG;

         /**
          * @brief Helicity tag name for CSCS file
          */
         static const std::string   HELICITYTAG;

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
