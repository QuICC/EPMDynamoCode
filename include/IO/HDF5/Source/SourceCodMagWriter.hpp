/** \file SourceCodMagWriter.hpp
 *  \brief Implementation of the source writer for the magnetic and codensity field 
 */

#ifndef SOURCECODMAGWRITER_HPP
#define SOURCECODMAGWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceCodWriter.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the source writer for the magnetic and codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodMagWriter: public SourceCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagWriter() {};
         
      protected:
         /**
          * @brief write magnetic source to file
          */
         virtual void writeMag();

         /**
          * @brief Pointer to the magnetic variable
          */
         const typename TSimTraits<TSimType>::MagType*  mpMagB;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> SourceCodMagWriter<TSimType, TSimTraits, TICSrc, TOCSrc>::SourceCodMagWriter(std::string name, SmartTruncation pTrunc)
      : SourceCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>(name, pTrunc), mpMagB(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodMagWriter<TSimType, TSimTraits, TICSrc, TOCSrc>::writeMag()
   {
      // Write magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->writeMagnetic(this->mpMagB->ocSrc().tor().data(), this->mpMagB->ocSrc().pol().data());
      }
   }

   /**
    * \brief Implementation of the source writer for the magnetic and codensity field. Specialisation for no OC source.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> class SourceCodMagWriter<TSimType, TSimTraits, TICSrc, false>: public SourceCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagWriter() {};
         
      protected:
         /**
          * @brief Write magnetic source from file
          */
         virtual void writeMag(){};

         /**
          * @brief Pointer to the magnetic variable
          */
         const typename TSimTraits<TSimType>::MagType*  mpMagB;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> SourceCodMagWriter<TSimType, TSimTraits, TICSrc, false>::SourceCodMagWriter(std::string name, SmartTruncation pTrunc)
      : SourceCodWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasCodICSource, TSimTraits<TSimType>::HasCodOCSource>(name, pTrunc)
   {
   }

}

#endif // SOURCECODMAGWRITER_HPP
