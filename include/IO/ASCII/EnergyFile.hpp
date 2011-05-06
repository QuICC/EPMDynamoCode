/** \file EnergyFile.hpp
 *  \brief Implementation of an erngy diagnostic file
 */

#ifndef ENERGYFILE_HPP
#define ENERGYFILE_HPP

// Configuration includes
//

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
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class EnergyFile: public ASCIIFieldWriter<TFieldType, ASCIIEWriter>
   {
      public:
         /**
         * @brief Constructor
         *
         * @param var Field variable
         * @param name File name
         * @param tsParams Timestep parameters
         */
         EnergyFile(TFieldType& var, std::string name, const TimestepParameters &tsParams, const int format = 0);

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

         /**
          * @brief Output format flag
          */
         const int mFormat;
   };

   template <typename TFieldType> EnergyFile<TFieldType>::EnergyFile(TFieldType &var, std::string name, const TimestepParameters &tsParams, int format)
      : ASCIIFieldWriter<TFieldType, ASCIIEWriter>(var, name + EnergyFileDefs::BASENAME, EnergyFileDefs::EXTENSION, EnergyFileDefs::HEADER, EnergyFileDefs::VERSION), mrTSParams(tsParams), mFormat(format)
   {
   }

   template <typename TFieldType> void EnergyFile<TFieldType>::write()
   {
      this->mrVar.rOc().updateSpectra();

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write energy
         this->mFile << this->mrTSParams.time() << "  " << this->mrVar.oc().energy().transpose();

         // If format = 1 output the energy harmonic degree spectrum evolution
         if(this->mFormat == 1)
         {
            for(int i = 1; i < this->mrVar.oc().spectrumL().cols(); ++i)
            {
               this->mFile << this->mrVar.oc().spectrumL().col(i).transpose();
            }
         }

         // Add newline at the end of ouput
         this->mFile << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

}

#endif // ENERGYFILE_HPP
