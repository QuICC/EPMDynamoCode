/** \file IndexControl.hpp
 */

#ifndef INDEXCONTROL_HPP
#define INDEXCONTROL_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/Typedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the operator index control
    */
   class IndexControl
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nOp  Number of operators in array
          */
         IndexControl(SmartSTrunc pSTrunc);

         /**
          * @brief Destructor
          */
         virtual ~IndexControl();

         /**
          * @brief Get the number of harmonic degrees
          */
         int nL() const;

         /**
          * @brief Get the total number of operators
          */
         int nOp() const;

         /**
          * @brief Set the total number of operators
          */
         int& rNOp();

         /**
          * @brief Get vector index for requested harmonic mode
          *
          * @param l Harmonic degree l
          */
         int opIdx(const int l) const;

         /**
          * @brief Set vector index for requested harmonic mode
          *
          * @param l Harmonic degree l
          */
         int& rOpIdx(const int l);

         /**
          * @brief Get array of different real orders
          */
         const ArrayI& realOrders(const int l) const;

         /**
          * @brief Set array of different real orders
          */
         ArrayI& rRealOrders(const int l);

         /**
          * @brief Get array of different imaginary orders
          */
         const ArrayI& imagOrders(const int l) const;

         /**
          * @brief Get array of different imaginary orders
          */
         ArrayI& rImagOrders(const int l);

         /**
          * @brief Has operator ?
          */
         bool hasOp(const int l) const;

         /**
          * @brief Has real harmonic orders ?
          */
         bool hasRealOrders(const int l) const;

         /**
          * @brief Has imaginary harmonic orders ?
          */
         bool hasImagOrders(const int l) const;

         /**
          * @brief Get array of harmonic degrees with defined operators
          */
         const ArrayI&  opLs() const;

         /**
          * @brief Create harmonic degrees list
          */
         void initLArray(bool makeAllOperators = false);
         
      protected:
      private:
         /**
          * @brief Number of harmonic degrees
          */
         int mNl;

         /**
          * @brief Number of operators
          */
         int mNop;

         /**
          * @brief Array of booleans to now if operator is present
          */
         ArrayB   mHasOp;

         /**
          * @brief Array of booleans to now if there are real orders
          */
         ArrayB   mHasRealOrders;

         /**
          * @brief Vector of booleans to now if imag orders are different
          */
         ArrayB   mHasImagOrders;

         /**
          * @brief Array of index for operators
          */
         ArrayI  mOpIdx;

         /**
          * @brief Array of harmonic degrees where operator is defined
          */
         ArrayI  mLArray;

         /**
          * @brief Array of required real harmonic orders
          */
         std::vector<ArrayI> mRealOrders;

         /**
          * @brief Array of required imaginary harmonic orders
          */
         std::vector<ArrayI> mImagOrders;

         /**
          * @brief Initialise the index arrays
          */
         void initOrdersInfos();
   };

   inline int IndexControl::nL() const
   {
      return this->mNl;
   }

   inline int IndexControl::nOp() const
   {
      return this->mNop;
   }

   inline int& IndexControl::rNOp()
   {
      return this->mNop;
   }

   inline int IndexControl::opIdx(const int l) const
   {
      return this->mOpIdx(l);
   }

   inline int& IndexControl::rOpIdx(const int l)
   {
      return this->mOpIdx(l);
   }

   inline const ArrayI& IndexControl::realOrders(const int l) const
   {
      return this->mRealOrders.at(l);
   }

   inline ArrayI& IndexControl::rRealOrders(const int l)
   {
      this->mHasRealOrders(l) = true;

      return this->mRealOrders.at(l);
   }

   inline const ArrayI& IndexControl::imagOrders(const int l) const
   {
      return this->mImagOrders.at(l);
   }

   inline ArrayI& IndexControl::rImagOrders(const int l)
   {
      this->mHasImagOrders(l) = true;

      return this->mImagOrders.at(l);
   }

   inline bool IndexControl::hasOp(const int l) const
   {
      return this->mHasOp(l);
   }

   inline bool IndexControl::hasRealOrders(const int l) const
   {
      return this->mHasRealOrders(l);
   }

   inline bool IndexControl::hasImagOrders(const int l) const
   {
      return this->mHasImagOrders(l);
   }

   inline const ArrayI& IndexControl::opLs() const
   {
      return this->mLArray;
   }

   /// typedef for a shared_ptr for the IndexControl
   typedef EPMDYNAMO_SHAREDPTRNS::shared_ptr<IndexControl>   SmartIdxControl;

}

#endif // INDEXCONTROL_HPP
