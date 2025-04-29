/** \file SpectralHDF5Reader.hpp
 *  \brief Implementation of a general spectral data file reader
 */

#ifndef SPECTRALHDF5READER_HPP
#define SPECTRALHDF5READER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/HDF5Reader.hpp"
#include "Domain/Truncation.hpp"
#include "Domain/SimulationTruncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the HDF5File class for the Statefile format
    */
   class SpectralHDF5Reader: public HDF5Reader
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
         SpectralHDF5Reader(std::string name, std::string ext, std::string header, std::string version, SmartTruncation   pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SpectralHDF5Reader() {};

         /**
          * @brief Write State to file
          */
         virtual void read() = 0;
         
      protected:
         
         /**
          * @brief Spectral truncatio information
          */
         SmartTruncation mpTrunc;

         /**
          * @brief Simulation truncation from file truncation
          */
         SmartSimTrunc   mpFileSim;

         /**
          * @brief Get run truncation information
          */
         SmartTruncation runTrunc() const;

         /**
          * @brief Get file truncation information
          */
         SmartSimTrunc fileTrunc() const;

         /**
          * @brief Set the read arguments
          */
         void setReadArguments();

         /**
          * @brief Read truncation information to file
          */
         void readTruncation();

         /**
          * @brief Check truncation compatibily between data and file
          */
         void checkTruncation();

      private:
         /**
          * @brief Set the maximum number of IO operations
          */
         void setMaxIO();
   };

   inline SmartTruncation   SpectralHDF5Reader::runTrunc()  const
   {
      return this->mpTrunc;
   }

   inline SmartSimTrunc   SpectralHDF5Reader::fileTrunc()  const
   {
      return this->mpFileSim;
   }

}

#endif // SPECTRALHDF5READER_HPP
