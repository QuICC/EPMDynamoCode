/** \file EnergyFile.hpp
 */

#ifndef ENERGYFILE_HPP
#define ENERGYFILE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IO/ASCII/EnergyFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of ASCIIFile for a general Energy ASCII file
    *
    * \param TSimType Type of the simulation
    * \param TFieldType Type of the field
    */
   template <typename TSimType, typename TFieldType> class EnergyFile: public ASCIIFieldWriter<TSimType, TFieldType, ASCIIEWriter>
   {
      public:
         /**
         * @brief Constructor
         *
         * @param var Field variable
         * @param name File name
         * @param tsParams Timestep parameters
         */
         EnergyFile(TFieldType& var, std::string name, const TimestepParameters &tsParams);

         /**
         * @brief Destructor
         */
         virtual ~EnergyFile() {};

         /**
          * @brief Write energy to file
          */
         virtual void write();
         
      protected:

      private:

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;
   };

   template <typename TSimType, typename TFieldType> EnergyFile<TSimType, TFieldType>::EnergyFile(TFieldType &var, std::string name, const TimestepParameters &tsParams)
      : ASCIIFieldWriter<TSimType, TFieldType, ASCIIEWriter>(var, name + EnergyFileDefs::BASENAME, EnergyFileDefs::EXTENSION, EnergyFileDefs::HEADER, EnergyFileDefs::VERSION), mrTSParams(tsParams)
   {
   }

   template <typename TSimType, typename TFieldType> void EnergyFile<TSimType, TFieldType>::write()
   {
      this->mrVar.rOc().updateSpectra();

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write energy
         this->mFile << this->mrTSParams.time() << "  " << this->mrVar.oc().energy().transpose() << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

}

#endif // ENERGYFILE_HPP
