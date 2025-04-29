/** \file SourceCodMagVelReader.hpp
 *  \brief Implementation of the source reader for the velocity and magnetic and codensity field 
 */

#ifndef SOURCECODMAGVELREADER_HPP
#define SOURCECODMAGVELREADER_HPP

// Configuration includes
//

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
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodMagVelReader: public SourceCodMagReader<TSimTraits, TSimTraits::HasMagICSource, TSimTraits::HasMagOCSource>
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
          * @brief Read velocity source from file
          */
         virtual void readVel();

         /**
          * @brief Pointer to the velocity variable
          */
         typename TSimTraits::VelType*  mpVelV;

      private:
   };

   template <typename TSimTraits, bool TICSrc, bool TOCSrc> SourceCodMagVelReader<TSimTraits, TICSrc, TOCSrc>::SourceCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : SourceCodMagReader<TSimTraits, TSimTraits::HasMagICSource, TSimTraits::HasMagOCSource>(name, pTrunc), mpVelV(NULL)
   {
   }

   template <typename TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodMagVelReader<TSimTraits, TICSrc, TOCSrc>::readVel()
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
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimTraits, bool TICSrc> class SourceCodMagVelReader<TSimTraits, TICSrc, false>: public SourceCodMagReader<TSimTraits, TSimTraits::HasMagICSource, TSimTraits::HasMagOCSource>
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
          * @brief Read magnetic source from file
          */
         virtual void readVel(){};

      private:
   };

   template <typename TSimTraits, bool TICSrc> SourceCodMagVelReader<TSimTraits, TICSrc, false>::SourceCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : SourceCodMagReader<TSimTraits, TSimTraits::HasMagICSource, TSimTraits::HasMagOCSource>(name, pTrunc)
   {
   }

}

#endif // SOURCECODMAGVELREADER_HPP
