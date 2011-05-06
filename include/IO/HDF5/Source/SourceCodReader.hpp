/** \file SourceCodReader.hpp
 *  \brief Implementation of the source reader for the codensity field 
 */

#ifndef SOURCECODREADER_HPP
#define SOURCECODREADER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceFileReaderBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the source reader for the codensity field 
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodReader: public SourceFileReaderBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodReader() {};
         
      protected:
         /**
          * @brief Read magnetic source from file
          */
         virtual void readCod();

         /**
          * @brief Pointer to the condensity variable
          */
         typename TSimTraits::CodType*  mpCodC;

      private:
   };

   template <typename TSimTraits, bool TICSrc, bool TOCSrc> SourceCodReader<TSimTraits, TICSrc, TOCSrc>::SourceCodReader(std::string name, SmartTruncation pTrunc)
      : SourceFileReaderBase(name, pTrunc), mpCodC(NULL)
   {
   }

   template <typename TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodReader<TSimTraits, TICSrc, TOCSrc>::readCod()
   {
      // Read codensity coefficients
      if(this->mpCodC != NULL)
      {
         this->readCodensity(this->mpCodC->rOcSrc().data());
      }
   }

   /**
    * \brief Implementation of the source reader for the codensity field. Specialisation for no OC source.
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimTraits, bool TICSrc> class SourceCodReader<TSimTraits, TICSrc, false>: public SourceFileReaderBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodReader() {};
         
      protected:
         /**
          * @brief Read magnetic source from file
          */
         virtual void readCod(){};

      private:
   };

   template <typename TSimTraits, bool TICSrc> SourceCodReader<TSimTraits, TICSrc, false>::SourceCodReader(std::string name, SmartTruncation pTrunc)
      : SourceFileReaderBase(name, pTrunc)
   {
   }

}

#endif // SOURCECODREADER_HPP
