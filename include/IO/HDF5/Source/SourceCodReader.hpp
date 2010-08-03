/** \file SourceCodReader.hpp
 *  \brief Implementation of the source reader for the codensity field 
 */

#ifndef SOURCECODREADER_HPP
#define SOURCECODREADER_HPP

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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodReader: public SourceFileReaderBase
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
         typename TSimTraits<TSimType>::CodType*  mpCodC;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> SourceCodReader<TSimType, TSimTraits, TICSrc, TOCSrc>::SourceCodReader(std::string name, SmartTruncation pTrunc)
      : SourceFileReaderBase(name, pTrunc), mpCodC(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodReader<TSimType, TSimTraits, TICSrc, TOCSrc>::readCod()
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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> class SourceCodReader<TSimType, TSimTraits, TICSrc, false>: public SourceFileReaderBase
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

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> SourceCodReader<TSimType, TSimTraits, TICSrc, false>::SourceCodReader(std::string name, SmartTruncation pTrunc)
      : SourceFileReaderBase(name, pTrunc)
   {
   }

}

#endif // SOURCECODREADER_HPP
