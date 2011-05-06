/** \file ImposedCodMagVelReader.hpp
 *  \brief Implementation of the imposed field reader for the velocity and magnetic and codensity field 
 */

#ifndef IMPOSEDCODMAGVELREADER_HPP
#define IMPOSEDCODMAGVELREADER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedCodMagReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the imposed field reader for the velocity and magnetic and codensity field 
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodMagVelReader: public ImposedCodMagReader<TSimTraits, TSimTraits::HasMagICImposed, TSimTraits::HasMagOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagVelReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagVelReader() {};
         
      protected:
         /**
          * @brief Read velocity imposed field from file
          */
         virtual void readVel();

         /**
          * @brief Pointer to the velocity variable
          */
         typename TSimTraits::VelType*  mpVelV;

      private:
   };

   template <typename TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodMagVelReader<TSimTraits, TICImposed, TOCImposed>::ImposedCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : ImposedCodMagReader<TSimTraits, TSimTraits::HasMagICImposed, TSimTraits::HasMagOCImposed>(name, pTrunc), mpVelV(NULL)
   {
   }

   template <typename TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodMagVelReader<TSimTraits, TICImposed, TOCImposed>::readVel()
   {
      // Read velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->readVelocity(this->mpVelV->rOc().rImposed().rTor().data(), this->mpVelV->rOc().rImposed().rPol().data());
      }
   }

   /**
    * \brief Implementation of the imposed field reader for the magnetic and codensity field. Specialisation for no OC imposed field.
    *
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    */
   template <typename TSimTraits, bool TICImposed> class ImposedCodMagVelReader<TSimTraits, TICImposed, false>: public ImposedCodMagReader<TSimTraits, TSimTraits::HasMagICImposed, TSimTraits::HasMagOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagVelReader(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagVelReader() {};
         
      protected:
         /**
          * @brief Read magnetic imposed field from file
          */
         virtual void readVel(){};

      private:
   };

   template <typename TSimTraits, bool TICImposed> ImposedCodMagVelReader<TSimTraits, TICImposed, false>::ImposedCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : ImposedCodMagReader<TSimTraits, TSimTraits::HasMagICImposed, TSimTraits::HasMagOCImposed>(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODMAGVELREADER_HPP
