/** \file ImposedCodMagReader.hpp
 *  \brief Implementation of the imposed field reader for the magnetic and codensity field 
 */

#ifndef IMPOSEDCODMAGREADER_HPP
#define IMPOSEDCODMAGREADER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedCodReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the imposed field reader for the magnetic and codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodMagReader: public ImposedCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagReader() {};
         
      protected:
         /**
          * @brief Read magnetic imposed field from file
          */
         virtual void readMag();

         /**
          * @brief Pointer to the magnetic variable
          */
         typename TSimTraits<TSimType>::MagType*  mpMagB;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodMagReader<TSimType, TSimTraits, TICImposed, TOCImposed>::ImposedCodMagReader(std::string name, SmartTruncation pTrunc)
      : ImposedCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>(name, pTrunc), mpMagB(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodMagReader<TSimType, TSimTraits, TICImposed, TOCImposed>::readMag()
   {
      // Read magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->readMagnetic(this->mpMagB->rOc().rImposed().rTor().data(), this->mpMagB->rOc().rImposed().rPol().data());
      }
   }

   /**
    * \brief Implementation of the imposed field reader for the magnetic and codensity field. Specialisation for no OC imposed field.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> class ImposedCodMagReader<TSimType, TSimTraits, TICImposed, false>: public ImposedCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagReader() {};
         
      protected:
         /**
          * @brief Read magnetic imposed field from file
          */
         virtual void readMag(){};

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> ImposedCodMagReader<TSimType, TSimTraits, TICImposed, false>::ImposedCodMagReader(std::string name, SmartTruncation pTrunc)
      : ImposedCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICImposed, TSimTraits<TSimType>::HasCodOCImposed>(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODMAGREADER_HPP
