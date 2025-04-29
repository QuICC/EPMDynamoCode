/** \file FileBase.hpp
 *  \brief Most basic part of a file implementation
 */

#ifndef FILEBASE_HPP
#define FILEBASE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief Most general features present in all file types
    */
   class FileBase
   {
      public:
         /**
         * @brief Constructor
         *
         * @param filename Name of the file
         * @param extension Extension of the file
         */
         FileBase(std::string filename, std::string extension);

         /**
         * @brief Destructor
         */
         virtual ~FileBase() {};

         /**
          * @brief Get filename
          */
         virtual std::string  filename() const;
         
      protected:

         /**
          * @brief Reset name
          *
          * @param name New name
          */
         void resetName(std::string name);

         /**
          * @brief Get the name
          */
         std::string  name() const;

         /**
          * @brief Get extension
          */
         std::string  extension() const;

         /**
          * @brief Does local core do IO?
          */
         bool doesIO() const;

         /**
          * @brief Rank of core
          */
         int rank() const;

         /**
          * @brief Rank allowed to do IO
          */
         int ioRank() const;

      private:

         /**
          * @brief Name of the file without extension
          */
         std::string mName;

         /**
          * @brief File extension
          */
         std::string mExt;

         /**
          * @brief Rank of the IO core
          */
         const int mcIORank;

         /**
          * @brief Rank of the local cpu
          */
         int mRank;

         /**
          * @brief Is this a parallel version of the code?
          */
         bool mIsParallel;

         /**
          * @brief Is local cpu doing IO?
          */
         bool mDoesIO;
   };

   inline std::string FileBase::filename() const
   {
      return this->mName + this->mExt;
   }

   inline std::string FileBase::name() const
   {
      return this->mName;
   }

   inline std::string FileBase::extension() const
   {
      return this->mExt;
   }

   inline bool FileBase::doesIO() const
   {
      return this->mDoesIO;
   }

   inline int FileBase::ioRank() const
   {
      return this->mcIORank;
   }

   inline int FileBase::rank() const
   {
      return this->mRank;
   }

}

#endif // FILEBASE_HPP
