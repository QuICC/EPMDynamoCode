/** \file ImposedCodMagVelWriter.hpp
 *  \brief Implementation of the imposed field writer for the velocity and magnetic and codensity field 
 */

#ifndef IMPOSEDCODMAGVELWRITER_HPP
#define IMPOSEDCODMAGVELWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedCodMagWriter.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the imposed field writer for the velocity and magnetic and codensity field 
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    * \tparam TOCImposed Has an OC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> class ImposedCodMagVelWriter: public ImposedCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagVelWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagVelWriter() {};
         
      protected:
         /**
          * @brief write velocity imposed field to file
          */
         virtual void writeVel();

         /**
          * @brief Pointer to the velocity variable
          */
         const typename TSimTraits<TSimType>::VelType*  mpVelV;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> ImposedCodMagVelWriter<TSimType, TSimTraits, TICImposed, TOCImposed>::ImposedCodMagVelWriter(std::string name, SmartTruncation pTrunc)
      : ImposedCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>(name, pTrunc), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed, bool TOCImposed> void ImposedCodMagVelWriter<TSimType, TSimTraits, TICImposed, TOCImposed>::writeVel()
   {
      // Write velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->writeVelocity(this->mpVelV->oc().imposed().tor().data(), this->mpVelV->oc().imposed().pol().data());
      }
   }

   /**
    * \brief Implementation of the imposed field writer for the velocity and magnetic and codensity field. Specialisation for no OC imposed field.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    * \tparam TICImposed Has an IC imposed field?
    */
   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> class ImposedCodMagVelWriter<TSimType, TSimTraits, TICImposed, false>: public ImposedCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>
   {
      public:
         
         /**
         * @brief Simple initialisation constructor (no work done)
         *
         * @param name Codensity variable
         * @param pTrunc Truncation information
         */
         ImposedCodMagVelWriter(std::string name, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedCodMagVelWriter() {};
         
      protected:
         /**
          * @brief Write velocity imposed field to file
          */
         virtual void writeVel(){};

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, bool TICImposed> ImposedCodMagVelWriter<TSimType, TSimTraits, TICImposed, false>::ImposedCodMagVelWriter(std::string name, SmartTruncation pTrunc)
      : ImposedCodMagWriter<TSimType, TSimTraits, TSimTraits<TSimType>::HasMagICImposed, TSimTraits<TSimType>::HasMagOCImposed>(name, pTrunc)
   {
   }

}

#endif // IMPOSEDCODMAGVELWRITER_HPP
