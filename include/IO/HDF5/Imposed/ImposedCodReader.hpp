/** \file ImposedCodReader.hpp
 *  \brief Implementation of the imposed field reader for the codensity field 
 */

#ifndef IMPOSEDCODREADER_HPP
#define IMPOSEDCODREADER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedFieldReaderBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the imposed field reader for the codensity field 
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodReader: public ImposedFieldReaderBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodReader() {};
         
      protected:
         /**
          * @brief Read magnetic imposed field from file
          */
         virtual void readCod();

         /**
          * @brief Pointer to the condensity variable
          */
         typename TSimTraits::CodType*  mpCodC;

      private:
   };

   template <typename TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodReader<TSimTraits, TICImposed, TOCImposed>::ImposedCodReader(std::string name, SmartTruncation pTrunc)
      : ImposedFieldReaderBase(name, pTrunc), mpCodC(NULL)
   {
   }

   template <typename TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodReader<TSimTraits, TICImposed, TOCImposed>::readCod()
   {
      // Read codensity coefficients
      if(this->mpCodC != NULL)
      {
         this->readCodensity(this->mpCodC->rOc().rImposed().data());
      }
   }

   /**
    * \brief Implementation of the imposed field reader for the codensity field. Specialisation for no OC imposed field.
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImpose Has an IC imposed field?
    */
   template <typename TSimTraits, bool TICImposed> class ImposedCodReader<TSimTraits, TICImposed, false>: public ImposedFieldReaderBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodReader() {};
         
      protected:
         /**
          * @brief Read magnetic imposed field from file
          */
         virtual void readCod(){};

      private:
   };

   template <typename TSimTraits, bool TICImposed> ImposedCodReader<TSimTraits, TICImposed, false>::ImposedCodReader(std::string name, SmartTruncation pTrunc)
      : ImposedFieldReaderBase(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODREADER_HPP
