/** \file ImposedCodMagWriter.hpp
 *  \brief Implementation of the imposed field writer for the magnetic and codensity field 
 */

#ifndef IMPOSEDCODMAGWRITER_HPP
#define IMPOSEDCODMAGWRITER_HPP

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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodMagWriter: public ImposedCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>
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
         const typename TSimTraits<TSimType>::MagType*  mpMagB;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodMagWriter<TSimType, TSimTraits, TICImposed, TOCImposed>::ImposedCodMagWriter(std::string name, SmartTruncation pTrunc)
      : ImposedCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>(name, pTrunc), mpMagB(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodMagWriter<TSimType, TSimTraits, TICImposed, TOCImposed>::writeMag()
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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> class ImposedCodMagWriter<TSimType, TSimTraits, TICImposed, false>: public ImposedCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>
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

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> ImposedCodMagWriter<TSimType, TSimTraits, TICImposed, false>::ImposedCodMagWriter(std::string name, SmartTruncation pTrunc)
      : ImposedCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODMAGWRITER_HPP
