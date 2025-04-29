/** \file ImposedCodMagWriter.hpp
 *  \brief Implementation of the imposed field writer for the magnetic and codensity field 
 */

#ifndef IMPOSEDCODMAGWRITER_HPP
#define IMPOSEDCODMAGWRITER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedCodWriter.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the imposed field writer for the magnetic and codensity field 
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodMagWriter: public ImposedCodWriter<TSimTraits, TSimTraits::HasCodICImposed, TSimTraits::HasCodOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagWriter() {};
         
      protected:
         /**
          * @brief write magnetic imposed field to file
          */
         virtual void writeMag();

         /**
          * @brief Pointer to the magnetic variable
          */
         const typename TSimTraits::MagType*  mpMagB;

      private:
   };

   template <typename TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodMagWriter<TSimTraits, TICImposed, TOCImposed>::ImposedCodMagWriter(std::string name, SmartTruncation pTrunc)
      : ImposedCodWriter<TSimTraits, TSimTraits::HasCodICImposed, TSimTraits::HasCodOCImposed>(name, pTrunc), mpMagB(NULL)
   {
   }

   template <typename TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodMagWriter<TSimTraits, TICImposed, TOCImposed>::writeMag()
   {
      // Write magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->writeMagnetic(this->mpMagB->oc().imposed().tor().data(), this->mpMagB->oc().imposed().pol().data());
      }
   }

   /**
    * \brief Implementation of the imposed field writer for the magnetic and codensity field. Specialisation for no OC impose field.
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    */
   template <typename TSimTraits, bool TICImposed> class ImposedCodMagWriter<TSimTraits, TICImposed, false>: public ImposedCodWriter<TSimTraits, TSimTraits::HasCodICImposed, TSimTraits::HasCodOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagWriter() {};
         
      protected:
         /**
          * @brief Write magnetic imposed field from file
          */
         virtual void writeMag(){};

      private:
   };

   template <typename TSimTraits, bool TICImposed> ImposedCodMagWriter<TSimTraits, TICImposed, false>::ImposedCodMagWriter(std::string name, SmartTruncation pTrunc)
      : ImposedCodWriter<TSimTraits, TSimTraits::HasCodICImposed, TSimTraits::HasCodOCImposed>(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODMAGWRITER_HPP
