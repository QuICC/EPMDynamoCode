/** \file IndexControl.cpp
 *  This file contains the definitions of the IndexControl class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "Operators/IndexControl.hpp"

// Project includes
//

namespace EPMDynamo {

   IndexControl::IndexControl(SmartSTrunc pSTrunc)
      : mNl(pSTrunc->local()->spec()->nL()), mNop(0), mHasOp(mNl), mHasRealOrders(mNl), mHasImagOrders(mNl), mOpIdx(mNl), mLArray()
   {
      // Initialise operator indexes to bad value
      this->mOpIdx.setConstant(-1);

      // Initialise operator presence check
      this->mHasOp.setConstant(false);

      // Initialise the harmonic orders information
      this->initOrdersInfos();
   }

   IndexControl::~IndexControl()
   {
   }

   void IndexControl::initLArray(bool makeAllOperators)
   {
      // Create temporary storage
      std::vector<int> tmp;

      // Initialise the harmonic orders indexes
      for(int l = 0; l < nL(); ++l)
      {
         if(this->mHasRealOrders(l) || this->mHasImagOrders(l))
         {
            tmp.push_back(l);
            this->mHasOp(l) = true;
         }
      }

      // Special case if all orders have to be present
      if(makeAllOperators)
      {
         this->mLArray.resize(nL());
         for(int l = 0; l < nL(); ++l)
         {
            this->mLArray(l) = l;
         }
      } else
      {
         this->mLArray.resize(tmp.size());
         for(int i = 0; i < this->mLArray.size(); ++i)
         {
            this->mLArray(i) = tmp.at(i);
         }
      }
   }

   void IndexControl::initOrdersInfos()
   {
      // Loop over the degrees
      for(int l = 0; l < nL(); ++l)
      {
         // Initialise bool to false
         this->mHasRealOrders(l) = false;
         this->mHasImagOrders(l) = false;

         // Create Array for storing indexes
         this->mRealOrders.push_back(ArrayI(1));
         this->mImagOrders.push_back(ArrayI(1));

         // Set it to a bad crashing value
         this->mRealOrders.at(l).setConstant(-1);
         this->mImagOrders.at(l).setConstant(-1);
      }
   }

}
