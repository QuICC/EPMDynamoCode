/** \file SourceWriterBase.hpp
 *  \brief Base of the source file writer
 */

#ifndef SOURCEWRITERBASE_HPP
#define SOURCEWRITERBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/SpectralHDF5NWriter.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the source file writer
    */
   class SourceWriterBase: public SpectralHDF5NWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name File name
         * @param pTrunc Truncation information
         */
         SourceWriterBase(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceWriterBase() {};

         /**
          * @brief Write State to file
          */
         virtual void write() = 0;

      protected:
         /**
          * @brief Create group for Codensity scalar and set mGroup
          *
          * @param scalar Codensity values
          */
         void writeCodensity(const std::vector<SpectralSHLShell> &scalar);

         /**
          * @brief Create group for Magnetic field and set mGroup
          *
          * @parm tor Toroidal magnetic values
          * @parm pol Poloidal magnetic values
          */
         void writeMagnetic(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol);

         /**
          * @brief Create group for Velocity field and set mGroup
          *
          * @parm tor Toroidal velocity values
          * @parm pol Poloidal velocity values
          */
         void writeVelocity(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol);

      private:
   };

}

#endif // SOURCEWRITERBASE_HPP
