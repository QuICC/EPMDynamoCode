/** \file HDF5NWriter.cpp
 *  \brief Source of the implementation of a numbering HDF5 file
 */

// System includes
//
#include <sstream>
#include <iomanip>

// External includes
//

// Class include
//
#include "IO/HDF5/HDF5NWriter.hpp"

// Project includes
//
#include "IO/HDF5/HDF5FileDefs.hpp"

namespace EPMDynamo {

   const int HDF5NWriter::msIDWidth = 4;

   HDF5NWriter::HDF5NWriter(std::string name, std::string ext, std::string header, std::string version)
      : HDF5Writer(name, ext, header, version), mCounter(0), mBaseName(name)
   {
   }

   void HDF5NWriter::updateName()
   {
      std::ostringstream   oss;

      oss << std::setfill('0') << std::setw(this->msIDWidth) << this->mCounter;

      this->resetName(this->mBaseName + oss.str());
   }

   void HDF5NWriter::init()
   {
   }

   void HDF5NWriter::preWrite()
   {
      // Update filename
      this->updateName();
      
      // Create new file
      this->open();
   }

   void HDF5NWriter::postWrite()
   {
      // Close the file
      this->close();

      // Increment counter
      ++this->mCounter;
   }

   void HDF5NWriter::finalise()
   {
   }

}
