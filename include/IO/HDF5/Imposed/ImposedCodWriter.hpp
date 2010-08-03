/** \file ImposedCodWriter.hpp
 *  \brief Implementation of the imposed field writer for the codensity field 
 */

#ifndef IMPOSEDCODWRITER_HPP
#define IMPOSEDCODWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedFieldWriterBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the imposed field writer for the codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodWriter: public ImposedFieldWriterBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodWriter() {};
         
      protected:
         /**
          * @brief Write codensity imposed field from file
          */
         virtual void writeCod();

         /**
          * @brief Pointer to the condensity variable
          */
         const typename TSimTraits<TSimType>::CodType*  mpCodC;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodWriter<TSimType, TSimTraits, TICImposed, TOCImposed>::ImposedCodWriter(std::string name, SmartTruncation pTrunc)
      : ImposedFieldWriterBase(name, pTrunc), mpCodC(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodWriter<TSimType, TSimTraits, TICImposed, TOCImposed>::writeCod()
   {
      // Write codensity coefficients
      if(this->mpCodC != NULL)
      {
         this->writeCodensity(this->mpCodC->oc().imposed().data());
      }
   }

   /**
    * \brief Implementation of the imposed field writer for the codensity field. Specialisation for no OC imposed field.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> class ImposedCodWriter<TSimType, TSimTraits, TICImposed, false>: public ImposedFieldWriterBase
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodWriter() {};
         
      protected:
         /**
          * @brief Write codensity imposed field from file
          */
         virtual void writeCod(){};

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> ImposedCodWriter<TSimType, TSimTraits, TICImposed, false>::ImposedCodWriter(std::string name, SmartTruncation pTrunc)
      : ImposedFieldWriterBase(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODWRITER_HPP
