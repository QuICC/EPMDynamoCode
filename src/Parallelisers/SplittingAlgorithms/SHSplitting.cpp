/** \file SHSplitting.cpp
 *  \brief Source of spherical harmonics load splitting algorithm
 */

// System includes
//
#include <utility>
#include <queue>
#include <set>
#include <boost/static_assert.hpp>

// External includes
//

// Class include
//
#include "Parallelisers/SplittingAlgorithms/SHSplitting.hpp"

// Project includes
//

namespace EPMDynamo {

   SHSplitting::SHSplitting(SmartSimTrunc pSimTrunc, const int nCore, const int id)
      : LoadSplitterBase(pSimTrunc, nCore, id)
   {
      BOOST_STATIC_ASSERT(false);
   }

   void SHSplitting::splitRTP(int &rR0, int &rNr, ArrayI &rTh0, ArrayI &rNth, const int id) const
   {
      this->splitR(rR0, rNr, id);

      rTh0.resize(rNr);
      rTh0.setConstant(0);

      int n = sim()->hoz()->nTh();
      rNth.resize(rNr);
      rNth.setConstant(n);
   }

   void SHSplitting::splitFDSH(ArrayI &rR0, ArrayI &rNr, ArrayI &rMs, const int id) const
   {
      rR0.resize(sim()->hoz()->nM());
      rR0.setConstant(0);

      int r0;
      int nR;
      this->splitR(r0, nR, id);

      rR0.resize(sim()->hoz()->nM());
      rR0.setConstant(r0);

      rNr.resize(sim()->hoz()->nM());
      rNr.setConstant(nR);

      rMs = sim()->hoz()->mArray();
   }

   void SHSplitting::splitSpec(ArrayI &rLs, std::vector<ArrayI> &rMs, const int id) const
   {
      int maxL;
      int maxM = this->sim()->hoz()->nM();
      if(this->sim()->hoz()->nM() % 2 == 0)
      {
         maxL = maxM;
      } else
      {
         maxL = maxM - 1;
      }

      std::queue<std::pair<int,int> >  modeQueue;
      for(int l=0; l < maxL/2; ++l)
      {
         for(int m=0; m < this->sim()->hoz()->nM(l); ++m)
         {
            modeQueue.push(std::make_pair(l, m));
         }

         for(int m=0; m < this->sim()->hoz()->nM(maxL-l-1); ++m)
         {
            modeQueue.push(std::make_pair(maxL-l-1, m));
         }
      }

      for(int l=maxL; l < this->sim()->hoz()->nL(); ++l)
      {
         for(int m=0; m < this->sim()->hoz()->nM(l); ++m)
         {
            modeQueue.push(std::make_pair(l, m));
         }
      }
      
      int h0, nH;
      balancedSplit(h0, nH, this->sim()->hoz()->nH(), id);

      for(int h=0; h < h0; ++h)
      {
         modeQueue.pop();
      }

      std::map<int, int>  lMap;
      std::multimap<int, int>  lMMap;
      for(int h=0; h < nH; ++h)
      {
         lMap.insert(modeQueue.front());
         lMMap.insert(modeQueue.front());
         modeQueue.pop();
      }

      std::map<int, int>::iterator  it;
      rLs.resize(lMap.size());
      int i = 0;
      rMs.clear();
      for(it=lMap.begin(); it != lMap.end(); it++, i++)
      {
         rLs(i) = (*it).first;
         rMs.push_back(ArrayI(lMMap.count(rLs(i))));
         for(int j = 0; j < rMs.at(i).size(); ++j)
         {
            rMs.at(i)(j) = (*it).second + j;
         }
      }

   }

}
