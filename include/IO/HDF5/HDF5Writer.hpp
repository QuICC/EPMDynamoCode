/** \file HDF5Writer.hpp
 *  \brief Implementation of a general HDF5 writer
 */

#ifndef HDF5WRITER_HPP
#define HDF5WRITER_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/HDF5TypesConv.hpp"
#include "IO/HDF5/HDF5File.hpp"

namespace EPMDynamo {

   /**
    * @brief Class responsible for writing data to an HDF5 file
    */
   class HDF5Writer: public HDF5File
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
         HDF5Writer(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~HDF5Writer() {};

         /**
          * @brief Initialise the file
          */
         virtual void init() = 0;

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalise() = 0;
         
      protected:
         /**
          * @brief Max collective IO Write operations over all CPUs
          */
         int mMaxCollIOWrite;

         /**
          * @brief Set total dataset size
          */
         virtual void setDatasetSize() = 0;

         /**
          * @brief Set dataset offsets
          */
         virtual void setDatasetOffsets() = 0;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Create the file info (Header and Version)
          */
         void createFileInfo();

         /**
          * @brief Write scalar dataset
          *
          * @param dsname HDF5 dataset name
          * @param val Scalar to write
          *
          * \tparam T Type of the scalar
          */
         template <typename T> void writeScalar(const std::string dsname, const T  val);

         /**
          * @brief Write Array dataset
          *
          * @param dsname HDF5 dataset name
          * @param inArray Array of values to write
          */
         void writeArray(const std::string dsname, const Array &inArray);

         /**
          * @brief Write vector of Matrix dataset
          *
          * @param dsname HDF5 dataset name
          * @param storage data to write
          */
         void writeMatrixVector(const std::string dsname, const std::vector<MatrixZ>  &storage);

         /**
          * @brief Write a 3D data vector dataset (for example the CSCS output)
          *
          * @param dsname HDF5 dataset name
          * @param storage Special data to write
          */
         void writeMatrixVector3D(const std::string dsname, const std::vector<Matrix>  &storage);

      private:
   };

   template <typename T> void HDF5Writer::writeScalar(const std::string dsname, const T val)
   {
      HDF5TypesConv genType;
      hid_t   type = genType.ptype<T>();

      hid_t dspace = H5Screate(H5S_SCALAR);

      mDataset = H5Dcreate(this->mGroup, dsname.c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      if(this->doesIO())
      {
         H5Dwrite(this->mDataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
      }

      H5Sclose(dspace);

      H5Dclose(this->mDataset);
   }

   /// Typedef for a smart pointer of a HDF5Writer
   typedef EPMSHARED_PTR<HDF5Writer> SmartHDF5Writer;

}

#endif // HDF5WRITER_HPP
