/** \file CSCSFileWriterBase.hpp
 *  \brief Base of the CSCS visualisation file format writer
 */

#ifndef CSCSFILEWRITERBASE_HPP
#define CSCSFILEWRITERBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/HDF5NWriter.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the CSCS visualisation file format writer
    */
   class CSCSFileWriterBase: public HDF5NWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param gridName Name of the grid group in HDF5 file
         * @param pTrunc Truncation information
         */
         CSCSFileWriterBase(const std::string gridName, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~CSCSFileWriterBase() {};

         /**
          * @brief Write State to file
          */
         virtual void write() = 0;

      protected:
         /**
          * @brief Spectral truncatio information
          */
         SmartTruncation mpTrunc;

         /**
          * @brief Get simulation truncation information
          */
         SmartTruncation runTrunc() const;

         /**
          * @brief Set the size of the dataset
          */
         void setDatasetSize();

         /**
          * @brief Set the offsets of the dataset
          */
         void setDatasetOffsets();

         /**
          * @brief Create the top data group (Grid name of CSCS file)
          */
         void createBaseGroup();

         /**
          * @brief Write Physical space grid to file
          *
          * @param radial Radial grid
          * @param theta Theta grid
          * @param phi Phi grid
          */
         void writeGrid(const Array& radial, const Array& theta, const Array& phi);

         /**
          * @brief Create group for Codensity scalar and set mGroup
          *
          * @param scalar Codensity values
          */
         void writeCodensity(const std::vector<SphericalShell> &scalar);

         /**
          * @brief Create group for Magnetic field and set mGroup
          *
          * @param radial Radial component of the magnetic field
          * @param theta Theta component of the magnetic field
          * @param phi Phi componenet of the magnetic field
          */
         void writeMagnetic(const std::vector<SphericalShell> &radial, const std::vector<SphericalShell> &theta, const std::vector<SphericalShell> &phi);

         /**
          * @brief Create group for Velocity field and set mGroup
          *
          * @param radial Radial component of the velocity field
          * @param theta Theta component of the velocity field
          * @param phi Phi componenet of the velocity field
          */
         void writeVelocity(const std::vector<SphericalShell> &radial, const std::vector<SphericalShell> &theta, const std::vector<SphericalShell> &phi);

      private:
         /**
          * @brief Name of the to group of the HDF5 file
          */
         const std::string    mGridName;

         /**
          * @brief Set the maximum number of IO operations
          */
         void setMaxIO();
   };

   inline SmartTruncation  CSCSFileWriterBase::runTrunc() const
   {
      return this->mpTrunc;
   }

}

#endif // CSCSFILEWRITERBASE_HPP
