/** \file HDF5Reader.hpp
 *  \brief Implementation of a general HDF5 file reader
 */

#ifndef HDF5READER_HPP
#define HDF5READER_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/HDF5TypesConv.hpp"
#include "IO/HDF5/HDF5File.hpp"

namespace EPMDynamo {

   /**
    * @brief Class responsible for reading data from HDF5 file
    */
   class HDF5Reader: public HDF5File
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name File name
         * @param ext File extension
         * @param header File header
         * @param version Version string
         */
         HDF5Reader(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~HDF5Reader() {};

         /**
          * @brief Initialise the file
          */
         void init();

         /**
          * @brief Read the content
          */

         virtual void read() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalise();
         
      protected:
         /**
          * @brief Number of blocks to read
          */
         hsize_t    mBlock;

         /**
          * @brief Max collective IO read operations over all CPUs
          */
         int mMaxCollIORead;

         /**
          * @brief Set data parameters
          */
         virtual void setReadArguments() = 0;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Check compatibility of opened file
          */
         void checkCompatibility();

         /**
          * @brief Read scalar dataset
          *
          * @param dsname HDF5 dataset name
          * @param val Storage for scalar to read
          *
          * \tparam T type of the scalar
          */
         template <typename T> void readScalar(const std::string dsname, T &val);

         /**
          * @brief Read vector of Matrix dataset
          *
          * @param dsname HDF5 dataset name
          * @param storage Storage for the data to read
          *
          * \tparam T type of the scalar
          */
         void readMatrixVector(const std::string dsname, std::vector<MatrixZ>  &storage);

      private:
   };

   template <typename T> void HDF5Reader::readScalar(const std::string dsname, T &val)
   {
      HDF5TypesConv genType;

      this->mDataset = H5Dopen(this->mGroup, dsname.c_str(), H5P_DEFAULT);
      H5Dread(this->mDataset, genType.ptype<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
      H5Dclose(this->mDataset);
   }

}

#endif // HDF5READER_HPP
