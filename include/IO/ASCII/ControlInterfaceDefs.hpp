/** \file ControlInterfaceDefs.hpp
 *  \brief Definitions and names used by the control interface
 */

#ifndef CONTROLINTERFACEDEFS_HPP
#define CONTROLINTERFACEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief Definitions for the user runtime interaction interface
    */
   class ControlInterfaceDefs
   {
      public:

         /**
          * @brief HEADER part for StdMessage file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for StdMessage file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of StdMessage file
          */
         static const std::string   NAME;

         /**
          * @brief EXTENSION of StdMessage file
          */
         static const std::string   EXTENSION;

         /**
          * @brief Run time input flag for stoping simulation
          */
         static const int STATUS_STOP;

         /**
          * @brief Run time input flag for keep going
          */
         static const int STATUS_GOON;

         /**
         * @brief Destructor
         */
         virtual ~ControlInterfaceDefs() {};
      private:
         /**
         * @brief Empty destructor
         */
         ControlInterfaceDefs();
   };



}

#endif // CONTROLINTERFACEDEFS_HPP
