/** \file SourceCodMagReader.hpp
 *  \brief Implementation of the source reader for the magnetic and codensity field 
 */

#ifndef SOURCECODMAGREADER_HPP
#define SOURCECODMAGREADER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceCodReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the source reader for the magnetic and codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodMagReader: public SourceCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagReader() {};
         
      protected:
         /**
          * @brief Read magnetic source from file
          */
         virtual void readMag();

         /**
          * @brief Pointer to the magnetic variable
          */
         typename TSimTraits<TSimType>::MagType*  mpMagB;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> SourceCodMagReader<TSimType, TSimTraits, TICSrc, TOCSrc>::SourceCodMagReader(std::string name, SmartTruncation pTrunc)
      : SourceCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>(name, pTrunc), mpMagB(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodMagReader<TSimType, TSimTraits, TICSrc, TOCSrc>::readMag()
   {
      // Read magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->readMagnetic(this->mpMagB->rOcSrc().rTor().data(), this->mpMagB->rOcSrc().rPol().data());
      }
   }

   /**
    * \brief Implementation of the source reader for the magnetic and codensity field. Specialisation for no OC source.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> class SourceCodMagReader<TSimType, TSimTraits, TICSrc, false>: public SourceCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagReader() {};
         
      protected:
         /**
          * @brief Read magnetic source from file
          */
         virtual void readMag(){};

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> SourceCodMagReader<TSimType, TSimTraits, TICSrc, false>::SourceCodMagReader(std::string name, SmartTruncation pTrunc)
      : SourceCodReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>(name, pTrunc)
   {
   }

}

#endif // SOURCECODMAGREADER_HPP
