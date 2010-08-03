/** \file ImposedCodMagVelReader.hpp
 *  \brief Implementation of the imposed field reader for the velocity and magnetic and codensity field 
 */

#ifndef IMPOSEDCODMAGVELREADER_HPP
#define IMPOSEDCODMAGVELREADER_HPP

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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodMagVelReader: public ImposedCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>
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
         typename TSimTraits<TSimType>::VelType*  mpVelV;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodMagVelReader<TSimType, TSimTraits, TICImposed, TOCImposed>::ImposedCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : ImposedCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>(name, pTrunc), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodMagVelReader<TSimType, TSimTraits, TICImposed, TOCImposed>::readVel()
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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> class ImposedCodMagVelReader<TSimType, TSimTraits, TICImposed, false>: public ImposedCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>
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

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> ImposedCodMagVelReader<TSimType, TSimTraits, TICImposed, false>::ImposedCodMagVelReader(std::string name, SmartTruncation pTrunc)
      : ImposedCodMagReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODMAGVELREADER_HPP
