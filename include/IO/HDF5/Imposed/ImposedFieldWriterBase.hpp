/** \file ImposedFieldWriterBase.hpp
 *  \brief Base of the HDF5 imposed field writers
 */

#ifndef IMPOSEDFIELDWRITERBASE_HPP
#define IMPOSEDFIELDWRITERBASE_HPP

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
    * \brief Base of the HDF5 imposed field writers
    */
   class ImposedFieldWriterBase: public SpectralHDF5NWriter
   {
      public:
         /**
         * @brief Constructor
         */
         ImposedFieldWriterBase(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedFieldWriterBase() {};

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
          * @param tor Toroidal magnetic values
          * @param pol Poloidal magnetic values
          */
         void writeMagnetic(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol);

         /**
          * @brief Create group for Velocity field and set mGroup
          *
          * @param tor Toroidal velocity values
          * @param pol Poloidal velocity values
          */
         void writeVelocity(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol);

      private:
   };

}

#endif // IMPOSEDFIELDWRITERBASE_HPP
