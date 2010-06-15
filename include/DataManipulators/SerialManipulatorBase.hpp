/** \file SerialManipulatorBase.hpp
 *  \brief Base for all the serial data manipulators.
 */

#ifndef SERIALMANIPULATORBASE_HPP
#define SERIALMANIPULATORBASE_HPP

// System includes
//
#include <queue>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "DataManipulators/ManipulatorBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Manipulation routines for  data
    *
    * \tparam TForward Type of data for the forward transform
    * \tparam TBackward Type of data for the backward transform
    */
   template <typename TForward, typename TBackward> class SerialManipulatorBase: public ManipulatorBase<TForward, TBackward>
   {
      public:
         /**
          * @brief Constructor
          *
          * \param pTrunc Pointer to the truncation
          * \param nFTmp Number of required forward storage units
          * \param nBTmp Number of required backward storage units
          */
         SerialManipulatorBase(const SmartTruncation pTrunc, const int nFTmp, const int nBTmp);

         /**
          * @brief Destructor
          */
         virtual ~SerialManipulatorBase() {};

         /**
          * @brief Release tempory storage after use
          * 
          * This routine doesn't imply putting it back, it just 
          * "releases" the control on it
          */
         void  releaseTmp(TForward &tmp);

         /**
          * @brief Release tempory storage after use
          * 
          * This routine doesn't imply putting it back, it just 
          * "releases" the control on it
          */
         void  releaseTmp(TBackward &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          */
         void  freeTmp(TForward &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          */
         void  freeTmp(TBackward &tmp);

         /**
          * @brief Send data
          *
          * \param data Data to send
          */
         void send(TForward &data);

         /**
          * @brief Send data
          *
          * \param data Data to send
          */
         void send(TBackward &data);

         /**
          * @brief Receive data
          *
          * \param rData Storage for received data
          */
         void receive(TForward &rData) = 0;

         /**
          * @brief Receive data
          *
          * \param rData Storage for received data
          */
         void receive(TBackward &rData) = 0;

         /**
          * @brief Set Desactivator state
          *
          * \param flag state
          */
         void setDesactivator(const bool flag);

         /**
          * @brief Finish initialisation
          */
         void setup();
         
      protected:
         /**
          * @brief Receive the sent data
          */
         TForward& receiveFData();

         /**
          * @brief Receive the sent data
          */
         TBackward& receiveBData();

      private:
         /**
          * @brief Temporary storage send queue for TForward
          */
         std::queue<TForward *> mFSendQueue;

         /**
          * @brief Temporary storage send queue for TBackward
          */
         std::queue<TBackward *> mBSendQueue;
   };

   template <typename TForward, typename TBackward> inline void SerialManipulatorBase<TForward, TBackward>::releaseTmp(TForward &tmp)
   {
   }

   template <typename TForward, typename TBackward> inline void SerialManipulatorBase<TForward, TBackward>::releaseTmp(TBackward &tmp)
   {
   }

   template <typename TForward, typename TBackward> inline void SerialManipulatorBase<TForward, TBackward>::setDesactivator(const bool flag)
   {
   }

   template <typename TForward, typename TBackward> inline void SerialManipulatorBase<TForward, TBackward>::freeTmp(TForward &tmp)
   {
      this->mFTmpQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void SerialManipulatorBase<TForward, TBackward>::freeTmp(TBackward &tmp)
   {
      this->mBTmpQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> SerialManipulatorBase<TForward, TBackward>::SerialManipulatorBase(const SmartTruncation pTrunc, const int nFTmp, const int nBTmp)
      : ManipulatorBase<TForward, TBackward>(pTrunc, nFTmp, nBTmp)
   {
   }

   template <typename TForward, typename TBackward> void SerialManipulatorBase<TForward, TBackward>::setup()
   {
   }

   template <typename TForward, typename TBackward> void SerialManipulatorBase<TForward, TBackward>::send(TForward &data)
   {
      this->mFSendQueue.push(&data);

      this->mPacksCounter--;
   }

   template <typename TForward, typename TBackward> void SerialManipulatorBase<TForward, TBackward>::send(TBackward &data)
   {
      this->mBSendQueue.push(&data);

      this->mPacksCounter--;
   }

   template <typename TForward, typename TBackward> TForward& SerialManipulatorBase<TForward, TBackward>::receiveFData()
   {
      this->mpFTmp = this->mFSendQueue.front();
      this->mFSendQueue.pop();

      return *this->mpFTmp;
   }

   template <typename TForward, typename TBackward> TBackward&  SerialManipulatorBase<TForward, TBackward>::receiveBData()
   {
      this->mpBTmp = this->mBSendQueue.front();
      this->mBSendQueue.pop();

      return *this->mpBTmp;
   }

}

#endif // SERIALMANIPULATORBASE_HPP
