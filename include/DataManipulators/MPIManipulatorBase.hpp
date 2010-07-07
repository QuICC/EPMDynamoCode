/** \file MPIManipulatorBase.hpp
 *  \brief Base for all the MPI data manipulators
 */

#ifndef MPIMANIPULATORBASE_HPP
#define MPIMANIPULATORBASE_HPP

// System includes
//
#include <map>
#include <set>
#include <vector>
#include <mpi.h>

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
    *
    * \epmBug Review/Cleanup code
    */
   template <typename TForward, typename TBackward> class MPIManipulatorBase: public ManipulatorBase<TForward, TBackward>
   {
      /// Typedef for the Forward transform data type
      typedef TForward  ForwardType;

      /// Typedef for the Backward transform data type
      typedef TBackward  BackwardType;

      public:
         /**
          * @brief Constructor
          *
          * \param pTrunc Pointer to the truncation
          * \param nFTmp Number of required forward storage units
          * \param nBTmp Number of required backward storage units
          */
         MPIManipulatorBase(const SmartTruncation pTrunc, const int nFTmp, const int nBTmp);

         /**
          * @brief Destructor
          */
         virtual ~MPIManipulatorBase();

         /**
          * @brief Release tempory storage after use
          * 
          * This routine doesn't imply putting it back, it just 
          * "releases" the control on it
          *
          * \param tmp Storage to release
          */
         void  releaseTmp(TForward &tmp);

         /**
          * @brief Release tempory storage after use
          * 
          * This routine doesn't imply putting it back, it just 
          * "releases" the control on it
          *
          * \param tmp Storage to release
          */
         void  releaseTmp(TBackward &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * \param tmp Storage to free (put back in queue)
          */
         void  freeTmp(TForward &tmp);

         /**
          * @brief Free tempory storage after use and put back into queue
          *
          * \param tmp Storage to free (put back in queue)
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
          *
          * \epmBug Removed minL(), might be bad idea
          */
         void receive(TForward &rData);

         /**
          * @brief Receive data
          *
          * \param rData Storage for received data
          *
          * \epmBug Removed minL(), might be bad idea
          */
         void receive(TBackward &rData);

         /**
          * @brief Start of a grouped "manual" Backward send
          *
          * \param entry Entry status for multiple entry
          */
         void initiateGroupedBSend(const int entry);

         /**
          * @brief Start of a grouped "manual" Forward send
          *
          * \param entry Entry status for multiple entry
          */
         void initiateGroupedFSend(const int entry);

         /**
          * @brief Set the interleaved stage entries ids
          *
          * \param fEntry Forward interleaved stage entry ID
          * \param bEntry Backward interleaved stage entry ID
          */
         void setInterStageIDs(const int fEntry, const int bEntry);

         /**
          * @brief Set the entry shift value
          *
          * \param shift Amount to shift the standard entry
          */
         void setEntryShift(const int shift);

         /**
          * @brief Set Desactivator for entries
          *
          * \param flag State of the desactivator
          */
         void setDesactivator(const bool flag);

         /**
          * @brief Finish initialisation
          */
         void setup();
         
      protected:
         /**
          * @brief Storage for the Forward datatypes
          */
         std::vector<std::map<TForward *, MPI_Datatype> >  mFTypes;

         /**
          * @brief Storage for the Backward datatypes
          */
         std::vector<std::map<TBackward *, MPI_Datatype> > mBTypes;

         /**
          * @brief Storage for the Forward buffers
          */
         std::vector<char *>  mFBuffers;

         /**
          * @brief Storage for the Backward buffers
          */
         std::vector<char *>  mBBuffers;

         /**
          * @brief Storage for the Forward buffer sizes
          */
         std::vector<int>  mFSizes;

         /**
          * @brief Storage for the Backward buffer sizes
          */
         std::vector<int>  mBSizes;

         /**
          * @brief Storage for the receive position pointers
          */
         std::vector<int>  mRecvPositions;

         /**
          * @brief Storage for the send position pointers
          */
         std::vector<int>  mSendPositions;

         /**
          * @brief Storage for the non blocking communication requests: Recv B
          */
         std::map<int, std::vector<MPI_Request> >  mRecvBRequests;

         /**
          * @brief Storage for the non blocking communication requests: Recv F
          */
         std::map<int, std::vector<MPI_Request> >  mRecvFRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send B
          */
         std::map<int, std::vector<MPI_Request> >  mSendBRequests;

         /**
          * @brief Storage for the non blocking communication requests: Send F
          */
         std::map<int, std::vector<MPI_Request> >  mSendFRequests;

         /**
          * @brief Storage for the index of the CPU group B
          */
         std::vector<int>  mCPUGroupB;

         /**
          * @brief Storage for the index of the CPU group F
          */
         std::vector<int>  mCPUGroupF;

         /**
          * @brief Get size of CPU group B
          */
         unsigned int sizeGroupB() const;

         /**
          * @brief Get size of CPU group B
          */
         unsigned int sizeGroupF() const;

         /**
          * @brief Get id mapping for CPU in group B
          */
         int idGroupB(const int id) const;

         /**
          * @brief Get id mapping for CPU in group B
          */
         int idGroupF(const int id) const;

         /**
          * @brief Get a pointer to the Recv B requests
          *
          * \param size "Size" of the requested request
          */
         MPI_Request * pRecvBRequests(const int size);

         /**
          * @brief Get a pointer to the Recv F requests
          *
          * \param size "Size" of the requested request
          */
         MPI_Request * pRecvFRequests(const int size);

         /**
          * @brief Get a pointer to the Send B requests
          *
          * \param size "Size" of the requested request
          */
         MPI_Request * pSendBRequests(const int size);

         /**
          * @brief Get a pointer to the Send F requests
          *
          * \param size "Size" of the requested request
          */
         MPI_Request * pSendFRequests(const int size);

         /**
          * @brief Boolean to trigger barrier for send to finish
          */
         bool  mIsSending;

         /**
          * @brief Boolean to trigger barrier for recv to finish
          */
         bool  mIsReceiving;

         /**
          * @brief Used value for desactivation
          */
         const int mDesactivationValue;

         /**
          * @brief Entry id for the interleaved stage for forward step
          */
         int mInterFEntry;

         /**
          * @brief Entry id for the interleaved stage for backward step
          */
         int mInterBEntry;

         /**
          * @brief Used to desactivate the entry
          */
         int mDesactivateEntry;

         /**
          * @brief Shift value for the entry
          */
         int mShiftEntry;

         /**
          * @brief Initialise the MPI datatypes
          */
         virtual void initTypes() = 0;

         /**
          * @brief Initialise the buffer sizes
          */
         void initBufferSizes();

         /**
          * @brief Initialise the positions
          */
         void initPositions();

         /**
          * @brief Setup the MPI communication requests requests
          */
         void setupRequests();

         /**
          * @brief Clear memory used by types
          */
         void cleanTypes();

         /**
          * @brief Clear memory used by buffers
          */
         void cleanBuffers();

         /**
          * @brief Clear memory used by requests
          */
         void cleanRequests();

         /**
          * @brief Build an TForward MPI Datatype
          */
         virtual void buildFType(TForward &data, MPI_Datatype &type, const int coreID) = 0;

         /**
          * @brief Build an TBackward MPI Datatype
          */
         virtual void buildBType(TBackward &data, MPI_Datatype &type, const int coreID) = 0;

         /**
          * @brief Does it need to wait to send the packed data?
          */
         bool isSending() const;

         /**
          * @brief Does it need to wait to receive the packed data?
          */
         bool isReceiving() const;

         /**
          * @brief Reset Receive positions
          */
         void resetRecvPositions();

         /**
          * @brief Reset Send positions
          */
         void resetSendPositions();

         /**
          * @brief Size of the forward packet
          */
         int sizeFPacket(const int id) const;

         /**
          * @brief Size of the backward packet
          */
         int sizeBPacket(const int id) const;

         /**
          * @brief Do send/receive loop for forward data
          */
         void initiateFSend();

         /**
          * @brief Do send/receive loop for backward data
          */
         void initiateBSend();

         /**
          * @brief Allocate the buffers
          */
         void allocateBuffers(const int maxPacks);

         /**
          * @brief Get ring recv source for id
          */
         int  recvSrc(const int id, const int ref, const int size) const;

         /**
          * @brief Get ring send destination for id
          */
         int  sendDest(const int id, const int ref, const int size) const;

      private:
   };

   template <typename TForward, typename TBackward> inline unsigned int MPIManipulatorBase<TForward, TBackward>::sizeGroupB() const
   {
      return this->mCPUGroupB.size();
   }

   template <typename TForward, typename TBackward> inline unsigned int MPIManipulatorBase<TForward, TBackward>::sizeGroupF() const
   {
      return this->mCPUGroupF.size();
   }

   template <typename TForward, typename TBackward> inline int MPIManipulatorBase<TForward, TBackward>::idGroupB(const int id) const
   {
      return this->mCPUGroupB.at(id);
   }

   template <typename TForward, typename TBackward> inline int MPIManipulatorBase<TForward, TBackward>::idGroupF(const int id) const
   {
      return this->mCPUGroupF.at(id);
   }

   template <typename TForward, typename TBackward> inline MPI_Request * MPIManipulatorBase<TForward, TBackward>::pRecvBRequests(const int size)
   {
      return &(this->mRecvBRequests[size].front());
   }

   template <typename TForward, typename TBackward> inline MPI_Request * MPIManipulatorBase<TForward, TBackward>::pRecvFRequests(const int size)
   {
      return &(this->mRecvFRequests[size].front());
   }

   template <typename TForward, typename TBackward> inline MPI_Request * MPIManipulatorBase<TForward, TBackward>::pSendBRequests(const int size)
   {
      return &(this->mSendBRequests[size].front());
   }

   template <typename TForward, typename TBackward> inline MPI_Request * MPIManipulatorBase<TForward, TBackward>::pSendFRequests(const int size)
   {
      return &(this->mSendFRequests[size].front());
   }

   template <typename TForward, typename TBackward> inline void MPIManipulatorBase<TForward, TBackward>::releaseTmp(TForward &tmp)
   {
      this->mFTmpQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void MPIManipulatorBase<TForward, TBackward>::releaseTmp(TBackward &tmp)
   {
      this->mBTmpQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void MPIManipulatorBase<TForward, TBackward>::freeTmp(TForward &tmp)
   {
      this->mFTmpQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline void MPIManipulatorBase<TForward, TBackward>::freeTmp(TBackward &tmp)
   {
      this->mBTmpQueue.push(&tmp);
   }

   template <typename TForward, typename TBackward> inline bool MPIManipulatorBase<TForward, TBackward>::isSending() const
   {
      return this->mIsSending;
   }

   template <typename TForward, typename TBackward> inline bool MPIManipulatorBase<TForward, TBackward>::isReceiving() const
   {
      return this->mIsReceiving;
   }

   template <typename TForward, typename TBackward> inline int MPIManipulatorBase<TForward, TBackward>::sizeFPacket(const int id) const
   {
      return this->mPacks*this->mFSizes.at(id);
   }

   template <typename TForward, typename TBackward> inline int MPIManipulatorBase<TForward, TBackward>::sizeBPacket(const int id) const
   {
      return this->mPacks*this->mBSizes.at(id);
   }

   template <typename TForward, typename TBackward> MPIManipulatorBase<TForward, TBackward>::MPIManipulatorBase(const SmartTruncation pTrunc, const int nFTmp, const int nBTmp)
      : ManipulatorBase<TForward, TBackward>(pTrunc, nFTmp, nBTmp), mIsSending(false), mIsReceiving(false), mDesactivationValue(-4242), mInterFEntry(-1), mInterBEntry(-1), mDesactivateEntry(mDesactivationValue), mShiftEntry(0)
   {
   }

   template <typename TForward, typename TBackward> MPIManipulatorBase<TForward, TBackward>::~MPIManipulatorBase()
   {
      // Cleanup Types memory
      this->cleanTypes();

      // Cleanup Buffers memory
      this->cleanBuffers();

      // Cleanup Requests memory
      this->cleanRequests();
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::setInterStageIDs(const int fEntry, const int bEntry)
   {
      this->mInterFEntry = fEntry;

      this->mInterBEntry = bEntry;
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::setEntryShift(const int shift)
   {
      this->mShiftEntry = shift;
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::setDesactivator(const bool flag)
   {
      if(flag)
      {
         this->mDesactivateEntry = this->mDesactivationValue;
      } else
      {
         this->mDesactivateEntry = 0;
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::cleanTypes()
   {
      // Cleanup the F Types
      typename std::vector<std::map<typename MPIManipulatorBase<TForward,TBackward>::ForwardType*, MPI_Datatype> >::iterator  itFT;
      typename std::map<typename MPIManipulatorBase<TForward,TBackward>::ForwardType*, MPI_Datatype>::iterator itF;

      for(itFT = this->mFTypes.begin(); itFT != this->mFTypes.begin(); ++itFT)
      {
         for(itF = (*itFT).begin(); itF != (*itFT).end(); ++itF)
         {
            MPI_Type_free(&((*itF).second));
         }
      }

      // Cleanup the B Types
      typename std::vector<std::map<typename MPIManipulatorBase<TForward,TBackward>::BackwardType*, MPI_Datatype> >::iterator  itBT;
      typename std::map<typename MPIManipulatorBase<TForward,TBackward>::BackwardType*, MPI_Datatype>::iterator  itB;

      for(itBT = this->mBTypes.begin(); itBT != this->mBTypes.begin(); ++itBT)
      {
         for(itB = (*itBT).begin(); itB != (*itBT).end(); ++itB)
         {
            MPI_Type_free(&((*itB).second));
         }
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::initBufferSizes()
   {
      int sze;
      int nCore = this->mpTrunc->para().nCore();

      for(int id = 0; id < nCore; ++id)
      {
         // Compute buffer sizes for F group
         MPI_Pack_size(1, this->mFTypes.at(id)[&(this->mFTmp.at(0))], MPI_COMM_WORLD, &sze);
         if(sze != 0)
         {
            this->mFSizes.push_back(sze);
            this->mCPUGroupF.push_back(id);
         }

         // Compute buffer sizes for F group
         MPI_Pack_size(1, this->mBTypes.at(id)[&(this->mBTmp.at(0))], MPI_COMM_WORLD, &sze);
         if(sze != 0)
         {
            this->mBSizes.push_back(sze);
            this->mCPUGroupB.push_back(id);
         }
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::allocateBuffers(const int maxPacks)
   {
      // Create B group buffers
      for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
      {
         this->mBBuffers.push_back(new char[maxPacks*this->mBSizes.at(id)]);
      }

      // Create F group buffers
      for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
      {
         this->mFBuffers.push_back(new char[maxPacks*this->mFSizes.at(id)]);
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::cleanBuffers()
   {
      // Freee memory of the B group buffers
      for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
      {
         delete[] this->mBBuffers.at(id);
      }

      // Freee memory of the F group buffers
      for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
      {
         delete[] this->mFBuffers.at(id);
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::cleanRequests()
   {
      std::map<int, std::vector<MPI_Request> >::iterator it;

      // Free requests from Recv B
      for(it = this->mRecvBRequests.begin(); it != this->mRecvBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Recv F
      for(it = this->mRecvFRequests.begin(); it != this->mRecvFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send B
      for(it = this->mSendBRequests.begin(); it != this->mSendBRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }

      // Free requests from Send F
      for(it = this->mSendFRequests.begin(); it != this->mSendFRequests.end(); it++)
      {
         for(unsigned int i = 0; i < (*it).second.size(); ++i)
         {
            MPI_Request_free(&((*it).second.at(i)));
         }
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::initPositions()
   {
      // Get maximum position size
      int maxSize = std::max(this->sizeGroupB(), this->sizeGroupF());

      // Initialise the position values
      for(int i = 0; i < maxSize; ++i)
      {
         this->mRecvPositions.push_back(0);

         this->mSendPositions.push_back(0);
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::setupRequests()
   {
      // Split behaviour between dynamic packs and collective communication
      int i0B;
      int i0F;
      if(this->isDynamicPacks())
      {
         i0B = 1;
         i0F = 1;
      } else
      {
         i0B = this->mMaxBPacks;
         i0F = this->mMaxFPacks;
      }

      std::vector<MPI_Request> tmp;

      int grpMe;
      int dest;
      int grpDest;
      int src;
      int grpSrc;
      int tag;

      // Initialise F requests
      for(int packs = i0F; packs <= this->mMaxFPacks ; ++packs)
      {
         // Initialise Recv B with empty requests
         this->mRecvBRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, tmp));

         for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
         {
            this->mRecvBRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Initialise Send F with empty requests
         this->mSendFRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, tmp));

         for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
         {
            this->mSendFRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Create Recv B requests
         for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
         {
            grpMe = (*std::find(this->mCPUGroupB.begin(), this->mCPUGroupB.end(), this->mpTrunc->para().id()));
            grpSrc = this->recvSrc(id, grpMe, this->sizeGroupB());
            src = this->idGroupB(grpSrc);
            tag = src;
            MPI_Recv_init(this->mBBuffers.at(grpSrc), packs*this->mBSizes.at(grpSrc), MPI_PACKED, src, tag, MPI_COMM_WORLD, &(this->mRecvBRequests[packs].at(grpSrc)));
         }

         // Create Send F requests
         for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
         {
            tag = this->mpTrunc->para().id();
            grpMe = (*std::find(this->mCPUGroupF.begin(), this->mCPUGroupF.end(), tag));
            grpDest = this->sendDest(id, grpMe, this->sizeGroupF());
            dest = this->idGroupF(grpDest);
            MPI_Send_init(this->mFBuffers.at(grpDest), packs*this->mFSizes.at(grpDest), MPI_PACKED, dest, tag, MPI_COMM_WORLD, &(this->mSendFRequests[packs].at(grpDest)));
         }
      }

      // Initialise B requests
      for(int packs = i0B; packs <= this->mMaxBPacks; ++packs)
      {
         // Initialise Recv F
         this->mRecvFRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, tmp));

         for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
         {
            this->mRecvFRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Initialise Send B
         this->mSendBRequests.insert(std::make_pair<int, std::vector<MPI_Request> >(packs, tmp));

         for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
         {
            this->mSendBRequests[packs].push_back(MPI_REQUEST_NULL);
         }

         // Create Recv F requests
         for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
         {
            grpMe = (*std::find(this->mCPUGroupF.begin(), this->mCPUGroupF.end(), this->mpTrunc->para().id()));
            grpSrc = this->recvSrc(id, grpMe, this->sizeGroupF());
            src = this->idGroupF(grpSrc);
            tag = src;
            MPI_Recv_init(this->mFBuffers.at(grpSrc), packs*this->mFSizes.at(grpSrc), MPI_PACKED, src, tag, MPI_COMM_WORLD, &(this->mRecvFRequests[packs].at(grpSrc)));
         }

         // Create Send B requests
         for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
         {
            tag = this->mpTrunc->para().id();
            grpMe = (*std::find(this->mCPUGroupB.begin(), this->mCPUGroupB.end(), tag));
            grpDest = this->sendDest(id, grpMe, this->sizeGroupB());
            dest = this->idGroupB(grpDest);
            MPI_Send_init(this->mBBuffers.at(grpDest), packs*this->mBSizes.at(grpDest), MPI_PACKED, dest, tag, MPI_COMM_WORLD, &(this->mSendBRequests[packs].at(grpDest)));
         }
      }
   }

   template <typename TForward, typename TBackward> int MPIManipulatorBase<TForward, TBackward>::sendDest(const int id, const int ref, const int size) const
   {
      // Create send ring
      return ((id + 1 + ref) % size);
   }

   template <typename TForward, typename TBackward> int MPIManipulatorBase<TForward, TBackward>::recvSrc(const int id, const int ref, const int size) const
   {
      // Create recv ring
      return ((size - 1 - id + ref) % size);
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::setup()
   {
      // initialise the positions
      this->initPositions();

      // Allocate memory for buffers
      this->allocateBuffers(std::max(this->mMaxFPacks, this->mMaxBPacks));

      // initialise the requests
      this->setupRequests();
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::resetRecvPositions()
   {
      std::vector<int>::iterator it;

      for(it = this->mRecvPositions.begin(); it != this->mRecvPositions.end(); ++it)
      {
         (*it) = 0;
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::resetSendPositions()
   {
      std::vector<int>::iterator it;

      for(it = this->mSendPositions.begin(); it != this->mSendPositions.end(); ++it)
      {
         (*it) = 0;
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::send(TForward &data)
   {
      if(this->isSending())
      {
         MPI_Waitall(this->sizeGroupF(), this->pSendFRequests(this->mPacks), MPI_STATUSES_IGNORE);

         this->mIsSending = false;
      }

      for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
      {
         MPI_Pack(MPI_BOTTOM, 1, this->mFTypes.at(this->idGroupF(id))[&data], this->mFBuffers.at(id), this->sizeFPacket(id), &(this->mSendPositions.at(id)), MPI_COMM_WORLD);
      }
      this->mPacksCounter--;

      if(this->atRecvStage())
      {
         this->initiateFSend();
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::initiateFSend()
   {
      // Prepost the receive calls
      MPI_Startall(this->sizeGroupB(), this->pRecvBRequests(this->mPacks));
      this->resetRecvPositions();
      this->mIsReceiving = true;

      // Post non blocking send calls 
      MPI_Startall(this->sizeGroupF(), this->pSendFRequests(this->mPacks));
      this->resetSendPositions();
      this->mIsSending = true;
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::send(TBackward &data)
   {
      if(this->isSending())
      {
         MPI_Waitall(this->sizeGroupB(), this->pSendBRequests(this->mPacks), MPI_STATUSES_IGNORE);

         this->mIsSending = false;
      }

      for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
      {
         MPI_Pack(MPI_BOTTOM, 1, this->mBTypes.at(this->idGroupB(id))[&data], this->mBBuffers.at(id), this->sizeBPacket(id), &(this->mSendPositions.at(id)), MPI_COMM_WORLD);
      }
      this->mPacksCounter--;

      if(this->atRecvStage())
      {
         this->initiateBSend();
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::initiateBSend()
   {
      // Prepost the receive calls
      MPI_Startall(this->sizeGroupF(), this->pRecvFRequests(this->mPacks));
      this->resetRecvPositions();
      this->mIsReceiving = true;

      // Post non blocking send calls 
      MPI_Startall(this->sizeGroupB(), this->pSendBRequests(this->mPacks));
      this->resetSendPositions();
      this->mIsSending = true;
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::receive(TForward &rData)
   {
      if(this->isReceiving())
      {
         int keepWaiting = this->sizeGroupF();
         int count = 0;
         ArrayI   idx(this->sizeGroupF());
         while(keepWaiting != 0)
         {
            MPI_Waitsome(this->sizeGroupF(), this->pRecvFRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);

            for(int id = 0; id < count; ++id)
            {
               MPI_Unpack(this->mFBuffers.at(idx(id)), this->sizeFPacket(idx(id)), &(this->mRecvPositions.at(idx(id))), MPI_BOTTOM, 1, this->mFTypes.at(this->idGroupF(idx(id)))[&rData], MPI_COMM_WORLD);
            }

            keepWaiting -= count;
         }

         this->mIsReceiving = false;
      } else
      {
         for(unsigned int id = 0; id < this->sizeGroupF(); ++id)
         {
            MPI_Unpack(this->mFBuffers.at(id), this->sizeFPacket(id), &(this->mRecvPositions.at(id)), MPI_BOTTOM, 1, this->mFTypes.at(this->idGroupF(id))[&rData], MPI_COMM_WORLD);
         }
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::receive(TBackward &rData)
   {
      if(this->isReceiving())
      {
         int keepWaiting = this->sizeGroupB();
         int count = 0;
         ArrayI   idx(this->sizeGroupB());
         while(keepWaiting != 0)
         {
            MPI_Waitsome(this->sizeGroupB(), this->pRecvBRequests(this->mPacks), &count, idx.data(), MPI_STATUSES_IGNORE);

            for(int id = 0; id < count; ++id)
            {
               MPI_Unpack(this->mBBuffers.at(idx(id)), this->sizeBPacket(idx(id)), &(this->mRecvPositions.at(idx(id))), MPI_BOTTOM, 1, this->mBTypes.at(this->idGroupB(idx(id)))[&rData], MPI_COMM_WORLD);
            }

            keepWaiting -= count;
         }

         this->mIsReceiving = false;
      } else
      {
         for(unsigned int id = 0; id < this->sizeGroupB(); ++id)
         {
            MPI_Unpack(this->mBBuffers.at(id), this->sizeBPacket(id), &(this->mRecvPositions.at(id)), MPI_BOTTOM, 1, this->mBTypes.at(this->idGroupB(id))[&rData], MPI_COMM_WORLD);
         }
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::initiateGroupedBSend(const int entry)
   {
      // This case controls the "simple" collective communication setup: sets first entry
      if(entry == 0 + this->mDesactivateEntry + this->mShiftEntry)
      {
         // Set first entry flag
         this->mFirstEntry = true;

         // Set additional pack to avoid automatic send
         this->mPacksCounter = 1;

         // Set the number of packs to send
         this->mPacks = this->mMaxBPacks;
      // This case controls the "simple" collective communication setup: sets second entry
      } else if(entry == 1 + this->mDesactivateEntry + this->mShiftEntry)
      {
         // Put system in recv stage
         this->mPacksCounter = 0;

         // initiate the send/recv loops
         this->initiateBSend();

         // prepare for second entry
         this->mFirstEntry = false;
      // This case controls the "simple nested" collective communication setup: sets first entry
      //    This is only called when a "nested collective" call is inside simple direct setup
      } else if(entry == 0 + this->mShiftEntry)
      {
         // Set first entry flag
         this->mFirstEntry = true;

         // Put system in recv stage
         this->mPacksCounter = 0;
      // This case controls the "nested" collective communication setup: sets interleaved entry
      } else if(entry == this->mInterBEntry)
      {
         // Put system in interlevead stage
         this->mPacksCounter = this->mDesactivationValue;
      }
   }

   template <typename TForward, typename TBackward> void MPIManipulatorBase<TForward, TBackward>::initiateGroupedFSend(const int entry)
   {
      // This case controls the "nested" collective communication setup: sets interleaved entry
      if(entry == this->mInterFEntry)
      {
         // unset first entry flag
         this->mFirstEntry = false;

         // Put system in special "nested" stage
         this->mPacksCounter = this->mDesactivationValue;
      }
      // This case controls the "simple" collective communication setup: sets first entry
      else if(entry == this->mInterFEntry + 1 + this->mDesactivateEntry)
      {
         // set first entry flag
         this->mFirstEntry = true;

         // Set additional pack to avoid automatic send
         this->mPacksCounter = 1;

         // Set the number of packs to send
         this->mPacks = this->mMaxFPacks;
      // This case controls the "simple" collective communication setup: sets second entry
      } else if(entry == this->mInterFEntry + 2 + this->mDesactivateEntry)
      {
         // Put system in recv stage
         this->mPacksCounter = 0;

         // initiate the send/recv loops
         this->initiateFSend();

         // prepare for second entry
         this->mFirstEntry = false;
      }
      // This case controls the "simple nested" collective communication setup: sets first entry
      //    This is only called when a "nested collective" call is inside simple direct setup
      else if(entry == this->mInterFEntry + 1)
      {
         // set first entry flag
         this->mFirstEntry = true;

         // Put system in recv stage
         this->mPacksCounter = 0;
      }
   }

}

#endif // MPIMANIPULATORBASE_HPP
