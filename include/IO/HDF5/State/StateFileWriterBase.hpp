/** \file StateFileWriterBase.hpp
 *  \brief Base for the HDF5 state file writer
 */

#ifndef STATEFILEWRITERBASE_HPP
#define STATEFILEWRITERBASE_HPP

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
    * \brief Base for the HDF5 state file writer
    */
   class StateFileWriterBase: public SpectralHDF5NWriter
   {
      public:
         /**
         * @brief Constructor
         */
         StateFileWriterBase(SmartTruncation   pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~StateFileWriterBase() {};

         /**
          * @brief Write State to file
          */
         virtual void write() = 0;

      protected:
         /**
          * @brief Write Physical parameters to file
          *
          * @param type Type of the parameters
          * @param params All physical parameters
          */
         void writePhysical(const std::string type, const Array& params);

         /**
          * @brief Write run information to file
          *
          * @param time Reached simulation time
          * @param step Current timestep
          */
         void writeRun(const EPMFloat time, const EPMFloat step);

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

#endif // STATEFILEWRITERBASE_HPP
