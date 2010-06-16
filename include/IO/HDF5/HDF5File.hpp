/** \file HDF5File.hpp
 *  \brief Implementation of a general HDF5 file
 */

#ifndef HDF5FILE_HPP
#define HDF5FILE_HPP

// System includes
//
#include <vector>
#include <hdf5.h>

// External includes
//

// Project includes
//
#include "IO/FileBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of FileBase for a general HDF5 file
    */
   class HDF5File: public FileBase
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name File name
         * @param ext File extension
         * @param header File header
         * @param version File version string
         * @param flags Special hdf5 flags
         */
         HDF5File(std::string name, std::string ext, std::string header, std::string version, const unsigned int flags);

         /**
         * @brief Destructor
         */
         virtual ~HDF5File() {};
         
      protected:
         /**
          * @brief Header of the file to check compatibility
          */
         std::string mHeader;

         /**
          * @brief Version of the file to check compatibility
          */
         std::string mVersion;

         /**
          * @brief File creation flag
          */
         unsigned int mFlags;

         /**
          * @brief Dataset property list
          */
         hid_t mDatasetPList;

         /**
          * @brief File creation property list
          */
         hid_t mFilePList;

         /**
          * @brief Handle to the actual HDF5 file
          */
         hid_t  mFile;

         /**
          * @brief HDF5 file group handle
          */
         hid_t mGroup;

         /**
          * @brief HDF5 file dataset handle
          */
         hid_t   mDataset;

         /**
          * @brief Size of the dimensions (HDF5 cols)
          */
         std::vector<hsize_t>  mFileDims;

         /**
          * @brief vector of offsets
          */
         std::vector< std::vector<hsize_t> > mFileOffsets;

         /**
          * @brief Set the dataset property list
          */
         void setDatasetPList();

         /**
          * @brief Free the dataset property list
          */
         void freeDatasetPList();

         /**
          * @brief Set the file property list
          */
         void setFilePList();

         /**
          * @brief Free the file property list
          */
         void freeFilePList();

      private:
   };

}

#endif // HDF5FILE_HPP
