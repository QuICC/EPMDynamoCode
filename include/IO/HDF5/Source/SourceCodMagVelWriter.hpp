/** \file SourceCodMagVelWriter.hpp
 *  \brief Implementation of the source writer for the velocity and magnetic and codensity field 
 */

#ifndef SOURCECODMAGVELWRITER_HPP
#define SOURCECODMAGVELWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceCodMagWriter.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the source writer for the velocity and magnetic and codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    * \tparam TOCSrc Has an OC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> class SourceCodMagVelWriter: public SourceCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagVelWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagVelWriter() {};
         
      protected:
         /**
          * @brief write velocity soure to file
          */
         virtual void writeVel();

         /**
          * @brief Pointer to the velocity variable
          */
         const typename TSimTraits<TSimType>::VelType*  mpVelV;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> SourceCodMagVelWriter<TSimType, TSimTraits, TICSrc, TOCSrc>::SourceCodMagVelWriter(std::string name, SmartTruncation pTrunc)
      : SourceCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>(name, pTrunc), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc, bool TOCSrc> void SourceCodMagVelWriter<TSimType, TSimTraits, TICSrc, TOCSrc>::writeVel()
   {
      // Write velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->writeVelocity(this->mpVelV->ocSrc().tor().data(), this->mpVelV->ocSrc().pol().data());
      }
   }

   /**
    * \brief Implementation of the source writer for the velocity and magnetic and codensity field. Specialisation for no OC source.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICSrc Has an IC source?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> class SourceCodMagVelWriter<TSimType, TSimTraits, TICSrc, false>: public SourceCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         SourceCodMagVelWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SourceCodMagVelWriter() {};
         
      protected:
         /**
          * @brief Write velocity source to file
          */
         virtual void writeVel(){};

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICSrc> SourceCodMagVelWriter<TSimType, TSimTraits, TICSrc, false>::SourceCodMagVelWriter(std::string name, SmartTruncation pTrunc)
      : SourceCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICSource, TSimTraits<TSimType>::HasMagOCSource>(name, pTrunc)
   {
   }

}

#endif // SOURCECODMAGVELWRITER_HPP
