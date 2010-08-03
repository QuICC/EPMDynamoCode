/** \file SourceCodMagVelReader.hpp
 *  \brief Implementation of the source reader for the velocity and magnetic and codensity field 
 */

#ifndef SOURCECODMAGVELREADER_HPP
#define SOURCECODMAGVELREADER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceCodMagReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the source reader for the velocity and magnetic and codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodMagVelReader: public SourceCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagVelReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagVelReader() {};
         
      protected:
         /**
          * @brief Read velocity soure from file
          */
         virtual void readVel();

         /**
          * @brief Pointer to the velocity variable
          */
         typename TSimTraits<TSimType>::VelType*  mpVelV;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> SourceCodMagVelReader<TSimType, TSimTraits, TICSrc, TOCSrc>::SourceCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : SourceCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>(name, pTrunc), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodMagVelReader<TSimType, TSimTraits, TICSrc, TOCSrc>::readVel()
   {
      // Read velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->readVelocity(this->mpVelV->rOcSrc().rTor().data(), this->mpVelV->rOcSrc().rPol().data());
      }
   }

   /**
    * \brief Implementation of the source reader for the magnetic and codensity field. Specialisation for no OC source.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> class SourceCodMagVelReader<TSimType, TSimTraits, TICSrc, false>: public SourceCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagVelReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagVelReader() {};
         
      protected:
         /**
          * @brief Read magnetic soure from file
          */
         virtual void readVel(){};

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> SourceCodMagVelReader<TSimType, TSimTraits, TICSrc, false>::SourceCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : SourceCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>(name, pTrunc)
   {
   }

}

#endif // SOURCECODMAGVELREADER_HPP
