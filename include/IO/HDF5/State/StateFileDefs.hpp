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
          * @brief Physical parameters part for State file
          */
         static const std::string   PHYSICALTAG;

         /**
          * @brief Ekman number tag for state file
          */
         static const std::string   EKMANTAG;

         /**
          * @brief Omega tag for state file
          */
         static const std::string   OMEGATAG;

         /**
          * @brief Viscosity tag for state file
          */
         static const std::string   NUTAG;

         /**
          * @brief Elsasser number tag for state file
          */
         static const std::string   ELSASSERTAG;

         /**
          * @brief Thermal Ekman number tag for state file
          */
         static const std::string   THERMALEKMANTAG;

         /**
          * @brief Magnetic Ekman number tag for state file
          */
         static const std::string   MAGNETICEKMANTAG;

         /**
          * @brief Magnetic Prandtl number tag for State file
          */
         static const std::string   MAGNETICPRANDTLTAG;

         /**
          * @brief Prandtl number tag for State file
          */
         static const std::string   PRANDTLTAG;

         /**
          * @brief Roberts number tag for State file
          */
         static const std::string   ROBERTSTAG;

         /**
          * @brief Rayleigh number tag for State file
          */
         static const std::string   RAYLEIGHTAG;

         /**
          * @brief Magnetic Rossby number tag for State file
          */
         static const std::string   ROSSBYTAG;

         /**
          * @brief Run parameters part for State file
          */
         static const std::string   RUNTAG;

         /**
          * @brief Time tag for State file
          */
         static const std::string   RUNTIMETAG;

         /**
          * @brief Timestep tag for State file
          */
         static const std::string   RUNSTEPTAG;

         /**
          * @brief Toroidal component tag for State file
          */
         static const std::string   TOROIDALTAG;

         /**
          * @brief Poloidal component tag for State file
          */
         static const std::string   POLOIDALTAG;

         /**
          * @brief Codensity field tag for State file
          */
         static const std::string   CODENSITYTAG;

         /**
          * @brief Magnetic field tag for State file
          */
         static const std::string   MAGNETICTAG;

         /**
          * @brief Velocity field tag for State file
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

         /**
          * @brief Filter id for using only symmetric component of field
          */
         static const int  SymmetricOnly = 4;

         /**
          * @brief Filter id for using only anti-symmetric component of field
          */
         static const int  AntisymmetricOnly = 5;

      private:
         /**
         * @brief Empty destructor
         */
         StateFileDefs();
   };
}

#endif // STATEFILEDEFS_HPP
