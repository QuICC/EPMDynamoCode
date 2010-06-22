/** \file ManipulatorBase.hpp
 *  \brief Base for all the data manipulators.
 */

#ifndef MANIPULATORBASE_HPP
#define MANIPULATORBASE_HPP

// System includes
//
#include <queue>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Basic Manipulation routines for templated data
    *
    * \tparam TForward Type of data for the forward transform
    * \tparam TBackward Type of data for the backward transform
    */
   template <typename TForward, typename TBackward> class ManipulatorBase
   {
      public:
         /**
          * @brief Constructor
          *
          * \param pTrunc Pointer to the truncation
          * \param nFTmp Number of required forward storage units
          * \param nBTmp Number of required backward storage units
          */
         ManipulatorBase(const SmartTruncation pTrunc, const int nFTmp, const int nBTmp);

         /**
          * @brief Destructor
          */
         virtual ~ManipulatorBase() {};

         /**
          * @brief Provide temporary storage for forward transform
          */
         TForward&  provideFTmp();

         /**
          * @brief Provide temporary storage for backward transform
          */
         TBackward&  provideBTmp();

         /**
          * @brief Release tempory storage after use
          * 
          * This routine doesn't imply putting it back, it just 
          * "releases" the control on it
          *
          * \param tmp Storage to release
          */
         virtual void  releaseTmp(TForward &tmp) = 0;

         /**
          * @brief Release tempory storage after use
          * 
          * This routine doesn't imply putting it back, it just 
          * "releases" the control on it
          *
          * \param tmp Storage to release
          */
         virtual void  releaseTmp(TBackward &tmp) = 0;

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * \param tmp Storage to free (put back in queue)
          */
         virtual void  freeTmp(TForward &tmp) = 0;

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * \param tmp Storage to free (put back in queue)
          */
         virtual void  freeTmp(TBackward &tmp) = 0;

         /**
          * @brief Send data
          *
          * \param data Data to send
          */
         virtual void send(TForward &data) = 0;

         /**
          * @brief Send data
          *
          * \param data Data to send
          */
         virtual void send(TBackward &data) = 0;

         /**
          * @brief Receive data
          *
          * \param rData Storage for received data
          */
         virtual void receive(TForward &rData) = 0;

         /**
          * @brief Receive data
          *
          * \param rData Storage for received data
          */
         virtual void receive(TBackward &rData) = 0;

         /**
          * @brief is transform at Send stage?
          */
         bool atSendStage() const;

         /**
          * @brief is transform at Receive stage?
          */
         bool atRecvStage() const;

         /**
          * @brief is transform at special stage?
          */
         bool atInterStage() const;

         /**
          * @brief Add n packs to next communication
          *
          * \param n Number of data packs to add
          */
         void addPacks(const int n);

         /**
          *  @brief Is first entry in function?
          */
         bool isFirstEntry() const;

         /**
          * @brief Set data packs control state
          *
          * \param state State of packs control
          */
         void setDynamicPacks(bool state);

         /**
          * @brief Set the max number of data packs
          *
          * \param maxFPacks Maximum number of forward packs
          * \param maxBPacks Maximum number of backward packs
          */
         void setMaxPacks(const int maxFPacks, const int maxBPacks);
         
      protected:
         /**
          * @brief Truncation information
          */
         SmartTruncation mpTrunc;

         /**
          * @brief Number of forward tempory storage scalars
          */
         const int mNFTmp;

         /**
          * @brief Number of backward tempory storage scalars
          */
         const int mNBTmp;

         /**
          * @brief Temporary storage queue for TForward
          */
         std::queue<TForward *> mFTmpQueue;

         /**
          * @brief Temporary storage queue for TBackward
          */
         std::queue<TBackward *> mBTmpQueue;

         /**
          * @brief Vector of TForward temporary storage
          */
         std::vector<TForward>  mFTmp;

         /**
          * @brief Vector of TBackward temporary storage
          */
         std::vector<TBackward>  mBTmp;

         /**
          * @brief TForward pointer
          */
         TForward   *mpFTmp;

         /**
          * @brief TBackward pointer
          */
         TBackward   *mpBTmp;

         /**
          * @brief Number of Packs to send
          */
         int mPacks;

         /**
          * @brief Counter for the packs to send
          */
         int mPacksCounter;

         /**
          * @brief bool for first entry check
          */
         bool mFirstEntry;

         /**
          * @brief bool for packs control state
          */
         bool mDynamicPacks;

         /**
          * @brief Maximum number of forward packs to send
          */
         int mMaxFPacks;

         /**
          * @brief Maximum number of backward Packs to send
          */
         int mMaxBPacks;

         /**
          *  @brief Is number of data packs set dynamicly?
          */
         bool isDynamicPacks() const;

      private:
         /**
          * @brief Initialise temporary storage
          */
         void initStorage();

         /**
          * @brief Initialise temporary storage queues (available storage)
          */
         void initQueues();
   };

   template <typename TForward, typename TBackward> inline bool ManipulatorBase<TForward, TBackward>::atSendStage() const
   {
      return (this->mPacksCounter > 0);
   }

   template <typename TForward, typename TBackward> inline bool ManipulatorBase<TForward, TBackward>::atRecvStage() const
   {
      return (this->mPacksCounter == 0);
   }

   template <typename TForward, typename TBackward> inline bool ManipulatorBase<TForward, TBackward>::atInterStage() const
   {
      return (this->mPacksCounter < 0);
   }

   template <typename TForward, typename TBackward> inline bool ManipulatorBase<TForward, TBackward>::isFirstEntry() const
   {
      return this->mFirstEntry;
   }

   template <typename TForward, typename TBackward> inline void ManipulatorBase<TForward, TBackward>::setDynamicPacks(bool state)
   {
      this->mDynamicPacks = state;
   }

   template <typename TForward, typename TBackward> inline bool ManipulatorBase<TForward, TBackward>::isDynamicPacks() const
   {
      return this->mDynamicPacks;
   }

   template <typename TForward, typename TBackward> inline void ManipulatorBase<TForward, TBackward>::addPacks(const int n)
   {
      // Increment data packs counter
      this->mPacksCounter += n;

      if(this->isDynamicPacks())
      {
         // Set number of packs
         this->mPacks = n;
      }
   }

   template <typename TForward, typename TBackward> ManipulatorBase<TForward, TBackward>::ManipulatorBase(const SmartTruncation pTrunc, const int nFTmp, const int nBTmp)
      : mpTrunc(pTrunc), mNFTmp(nFTmp), mNBTmp(nBTmp), mPacks(0), mPacksCounter(0), mFirstEntry(true), mDynamicPacks(true), mMaxFPacks(0), mMaxBPacks(0)
   {
      // Init the temporary storage
      this->initStorage();

      // Init the temporary storage queues
      this->initQueues();
   }

   template <typename TForward, typename TBackward> void ManipulatorBase<TForward, TBackward>::setMaxPacks(const int maxFPacks, const int maxBPacks)
   {
      this->mMaxFPacks = maxFPacks;
      this->mMaxBPacks = maxBPacks;
   }

   template <typename TForward, typename TBackward> void ManipulatorBase<TForward, TBackward>::initStorage()
   {
      for(int i=0; i < this->mNFTmp; ++i)
      {
         this->mFTmp.push_back(TForward(this->mpTrunc));
      }

      for(int i=0; i < this->mNBTmp; ++i)
      {
         this->mBTmp.push_back(TBackward(this->mpTrunc));
      }
   }

   template <typename TForward, typename TBackward> void ManipulatorBase<TForward, TBackward>::initQueues()
   {
      for(int i=0; i < this->mNFTmp; ++i)
      {
         this->mFTmpQueue.push(&(this->mFTmp.at(i)));
      }

      for(int i=0; i < this->mNBTmp; ++i)
      {
         this->mBTmpQueue.push(&(this->mBTmp.at(i)));
      }
   }

   template <typename TForward, typename TBackward> TForward&  ManipulatorBase<TForward, TBackward>::provideFTmp()
   {
      this->mpFTmp = this->mFTmpQueue.front();
      this->mFTmpQueue.pop();

      return *this->mpFTmp;
   }

   template <typename TForward, typename TBackward> TBackward&  ManipulatorBase<TForward, TBackward>::provideBTmp()
   {
      this->mpBTmp = this->mBTmpQueue.front();
      this->mBTmpQueue.pop();

      return *this->mpBTmp;
   }

}

#endif // MANIPULATORBASE_HPP
