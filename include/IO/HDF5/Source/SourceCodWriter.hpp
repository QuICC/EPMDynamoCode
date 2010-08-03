/** \file SourceCodWriter.hpp
 *  \brief Implementation of the source writer for the codensity field 
 */

#ifndef SOURCECODWRITER_HPP
#define SOURCECODWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceFileWriterBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the source writer for the codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodWriter: public SourceFileWriterBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodWriter() {};
         
      protected:
         /**
          * @brief Write codensity source from file
          */
         virtual void writeCod();

         /**
          * @brief Pointer to the condensity variable
          */
         const typename TSimTraits<TSimType>::CodType*  mpCodC;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> SourceCodWriter<TSimType, TSimTraits, TICSrc, TOCSrc>::SourceCodWriter(std::string name, SmartTruncation pTrunc)
      : SourceFileWriterBase(name, pTrunc), mpCodC(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodWriter<TSimType, TSimTraits, TICSrc, TOCSrc>::writeCod()
   {
      // Write codensity coefficients
      if(this->mpCodC != NULL)
      {
         this->writeCodensity(this->mpCodC->ocSrc().data());
      }
   }

   /**
    * \brief Implementation of the source writer for the codensity field. Specialisation for no OC source.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> class SourceCodWriter<TSimType, TSimTraits, TICSrc, false>: public SourceFileWriterBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodWriter() {};
         
      protected:
         /**
          * @brief Write codensity source from file
          */
         virtual void writeCod(){};

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> SourceCodWriter<TSimType, TSimTraits, TICSrc, false>::SourceCodWriter(std::string name, SmartTruncation pTrunc)
      : SourceFileWriterBase(name, pTrunc)
   {
   }

}

#endif // SOURCECODWRITER_HPP
