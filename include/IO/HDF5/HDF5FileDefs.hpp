/** \file HDF5FileDefs.hpp
 *  \brief General defininitions and names for the HDF5 files
 */

#ifndef HDF5FILEDEFS_HPP
#define HDF5FILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief General defininitions and names for the HDF5 files
    */
   class HDF5FileDefs
   {
      public:

         /**
          * @brief HEADER part for State file
          */
         static const std::string   HEADERTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   VERSIONTAG;

         /**
         * @brief Destructor
         */
         virtual ~HDF5FileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         HDF5FileDefs();
   };
}

#endif // HDF5FILEDEFS_HPP
