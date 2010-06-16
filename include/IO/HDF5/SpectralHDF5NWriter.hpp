/** \file SpectralHDF5NWriter.hpp
 *  \brief Implementation of a general spectral data file writer
 */

#ifndef SPECTRALHDF5WRITER_HPP
#define SPECTRALHDF5WRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/Typedefs.hpp"
#include "IO/HDF5/HDF5NWriter.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the HDF5File class for a general spectral data format
    */
   class SpectralHDF5NWriter: public HDF5NWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name File name
         * @param ext File extension
         * @param header File header
         * @param version Version string
         * @param pTrunc Truncation information
         */
         SpectralHDF5NWriter(std::string name, std::string ext, std::string header, std::string version, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SpectralHDF5NWriter() {};

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
          * @brief Get run truncation information
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
          * @brief Write truncation information to file
          *
          * @param N Simulation wide radial truncation
          * @param L Simulation wide harmonic degree truncation
          * @param M Simulation wide harmonic order truncation
          * @param Mp Assumed longitudinal symmetry
          */
         void writeTruncation(const int N, const int L, const int M, const int Mp);

      private:
         /**
          * @brief Set the maximum number of IO operations
          */
         void setMaxIO();
   };

   inline SmartTruncation   SpectralHDF5NWriter::runTrunc() const
   {
      return this->mpTrunc;
   }

}

#endif // SPECTRALHDF5WRITER_HPP
