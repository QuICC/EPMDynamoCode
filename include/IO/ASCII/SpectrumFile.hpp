/** \file SpectrumFile.hpp
 *  \brief Implementation of an ASCII spectrum diagnostic file
 */

#ifndef SPECTRUMFILE_HPP
#define SPECTRUMFILE_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IO/ASCII/SpectrumFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIRWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of ASCIIFile for a general spectrum ASCII file
    *
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class SpectrumFile: public ASCIIFieldWriter<TFieldType, ASCIIRWriter>
   {
      public:
         /**
         * @brief Constructor
         *
         * @param var Field variable
         * @param name File name
         */
         SpectrumFile(TFieldType& var, std::string name);

         /**
         * @brief Destructor
         */
         virtual ~SpectrumFile() {};

         /**
          * @brief Write spectrum to file
          */
         virtual void write();
         
      protected:

      private:
   };

   template <typename TFieldType> SpectrumFile<TFieldType>::SpectrumFile(TFieldType &var, std::string name)
      : ASCIIFieldWriter<TFieldType, ASCIIRWriter>(var, name + SpectrumFileDefs::BASENAME, SpectrumFileDefs::EXTENSION, SpectrumFileDefs::HEADER, SpectrumFileDefs::VERSION)
   {
   }

   template <typename TFieldType> void SpectrumFile<TFieldType>::write()
   {
      this->mrVar.rOc().updateSpectra();

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write M spectrum
         this->mFile << "#M Spectrum" << std::endl;
         this->mFile << this->mrVar.oc().spectrumM() << std::endl;

         // Add newlines
         this->mFile << std::endl << std::endl;

         // Write L spectrum
         this->mFile << "#L Spectrum" << std::endl;
         this->mFile << this->mrVar.oc().spectrumL() << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

}

#endif // SPECTRUMFILE_HPP
