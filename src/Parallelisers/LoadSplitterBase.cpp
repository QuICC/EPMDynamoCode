/** \file LoadSplitterBase.cpp
 *  \brief Source of the base of the implementation of the load splitting algorithms
 */

// System includes
//
#include <algorithm>
#include <set>

// External includes
//

// Class include
//
#include "Parallelisers/LoadSplitterBase.hpp"

// Project includes
//
#include "General/EPMException.hpp"

namespace EPMDynamo {

   LoadSplitterBase::LoadSplitterBase(SmartSimTrunc pSimTrunc, const int nCore, const int id)
      : mpSimTrunc(pSimTrunc), mNcore(nCore), mID(id)
   {
   }

   void LoadSplitterBase::balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id) const
   {
      // Avoid splitting with zero elements
      if(tot < parts)
      {
         throw EPMException("LoadSplitterBase::balancedSplit", "Number of parts is bigger than total!");
      }

      if(parts != 1)
      {
         nN = 0;
         n0 = 0;
         for(int i = 0; i < tot; i++)
         {
            if(i % parts == id)
            {
               nN++;
            }
            else if(i % parts < id)
            {
               n0++;
            }
         }
      } else
      {
         n0 = 0;
         nN = tot;
      }
   }

   void LoadSplitterBase::balancedSplit(int &n0, int &nN, const int tot, const int id) const
   {
      // Avoid splitting with zero elements
      if(tot < this->nCore())
      {
         throw EPMException("LoadSplitterBase::balancedSplit", "Number of CPUs is bigger than total!");
      }

      nN = 0;
      n0 = 0;
      for(int i = 0; i < tot; i++)
      {
         if(i % this->nCore() == id)
         {
            nN++;
         }
         else if(i % this->nCore() < id)
         {
            n0++;
         }
      }
   }

   double LoadSplitterBase::computeRTPInbalance(SmartLocalTrunc pLocal) const
   {
      // Compute optimal load splitting
      double optimalLoad = static_cast<double>(this->sim()->rad()->nR()*this->sim()->hoz()->nTh())/static_cast<double>(this->nCore());

      // Compute implemented load
      double cpuLoad = 0.0;
      for(int n = 0; n < pLocal->rtp()->nR(); ++n)
      {
         cpuLoad += static_cast<double>(pLocal->rtp()->nTh(n));
      }

      // Get load inbalance
      double inbalance;
      inbalance = std::abs((cpuLoad - optimalLoad)/optimalLoad);

      return  inbalance;
   }

   double LoadSplitterBase::computeFDSHInbalance(SmartLocalTrunc pLocal) const
   {
      ArrayI   nLs = this->sim()->hoz()->nLArray();

      // Compute optimal load splitting
      double optimalLoad = 0.0;
      for(int n = 0; n < this->sim()->rad()->nR(); ++n)
      {
         for(int m =0; m < this->sim()->hoz()->nM(); m++)
         {
            optimalLoad += static_cast<double>(nLs(m));
         }
      }
      optimalLoad /= static_cast<double>(this->nCore()),


      // Compute implemented load
      nLs = pLocal->fdsh()->nLArray(this->sim()->hoz()->nL());

      double cpuLoad = 0.0;
      for(int m =0; m < pLocal->fdsh()->nM(); m++)
      {
         for(int n = 0; n < pLocal->fdsh()->nR(m); ++n)
         {
            cpuLoad += static_cast<double>(nLs(m));
         }
      }

      // Get load inbalance
      double inbalance;
      inbalance = std::abs((cpuLoad - optimalLoad)/optimalLoad);

      return inbalance;
   }

   double LoadSplitterBase::computeSpecInbalance(SmartLocalTrunc pLocal) const
   {
      // Compute optimal load splitting
      double optimalLoad = 0.0;
      for(int l = 0; l < this->sim()->hoz()->nL(); ++l)
      {
         optimalLoad += static_cast<double>(this->sim()->hoz()->nM(l));
      }
      optimalLoad /= static_cast<double>(this->nCore());

      // Compute implemented load
      double cpuLoad = 0.0;
      for(int l = 0; l < pLocal->spec()->nL(); ++l)
      {
         cpuLoad += static_cast<double>(pLocal->spec()->nM(l));
      }

      // Get load inbalance
      double inbalance;
      inbalance = std::abs((cpuLoad - optimalLoad)/optimalLoad);

      return inbalance;
   }

   bool LoadSplitterBase::checkRTPConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const
   {
      // Multiply indexes to get unique IDs
      double radMult = 1.0e6;
      double thMult = 1.0;

      // Create the set of IDs from simulation wide truncations
      std::set<double>   simSet;
      for(int n=0; n < this->sim()->rad()->nR(); ++n)
      {
         for(int th=0; th < this->sim()->hoz()->nTh(); ++th)
         {
            simSet.insert(n*radMult + th*thMult);
         }
      }

      // Create the set of IDs from simulation splitted truncations
      std::set<double>   remSet;
      int n_;
      int th_;
      for(int id=0; id < this->nCore(); ++id)
      {
         for(int n=0; n < pRemote.at(id)->rtp()->nR(); ++n)
         {
            n_ = n + pRemote.at(id)->rtp()->r0();
            for(int th=0; th < pRemote.at(id)->rtp()->nTh(n); ++th)
            {
               th_ = th + pRemote.at(id)->rtp()->th0(n);
               remSet.insert(n_*radMult + th_*thMult);
            }
         }
      }

      // Check if both sets contain the same data
      if(std::equal(simSet.begin(), simSet.end(), remSet.begin()))
      {
         return true;
      } else
      {
         return false;
      }
   }

   bool LoadSplitterBase::checkFDSHConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const
   {
      // Multiply indexes to get unique IDs
      double radMult = 1.0e6;
      double mMult = 1.0;

      // Create the set of IDs from simulation wide truncations
      std::set<double>   simSet;
      for(int n = 0; n < this->sim()->rad()->nR(); ++n)
      {
         for(int m =0; m < this->sim()->hoz()->nM(); m++)
         {
            simSet.insert(n*radMult + m*mMult);
         }
      }

      // Create the set of IDs from simulation splitted truncations
      std::set<double>   remSet;
      int n_;
      int m_;
      ArrayI   ms;
      for(int id=0; id < this->nCore(); ++id)
      {
         ms = pRemote.at(id)->fdsh()->mArray();
         for(int m = 0; m < pRemote.at(id)->fdsh()->nM(); m++)
         {
            m_ = ms(m);
            for(int n = 0; n < pRemote.at(id)->fdsh()->nR(m); ++n)
            {
               n_ = n + pRemote.at(id)->fdsh()->r0(m);
               remSet.insert(n_*radMult + m_*mMult);
            }
         }
      }

      // Check if both sets contain the same data
      if(std::equal(simSet.begin(), simSet.end(), remSet.begin()))
      {
         return true;
      } else
      {
         return false;
      }
   }

   bool LoadSplitterBase::checkSpecConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const
   {
      // Multiply indexes to get unique IDs
      double lMult = 1.0e6;
      double mMult = 1.0;

      // Create the set of IDs from simulation wide truncations
      std::set<double>   simSet;
      for(int l = 0; l < this->sim()->hoz()->nL(); ++l)
      {
         for(int m = 0; m < this->sim()->hoz()->nM(l); ++m)
         {
            simSet.insert(l*lMult + m*mMult);
         }
      }

      // Create the set of IDs from simulation splitted truncations
      std::set<double>   remSet;
      int l_;
      int m_;
      ArrayI   ls;
      for(int id=0; id < this->nCore(); ++id)
      {
         ls = pRemote.at(id)->spec()->lArray();
         for(int l = 0; l < pRemote.at(id)->spec()->nL(); l++)
         {
            l_ = ls(l);
            for(int m = 0; m < pRemote.at(id)->spec()->nM(l); ++m)
            {
               m_ = pRemote.at(id)->spec()->m(m,l);
               remSet.insert(l_*lMult + m_*mMult);
            }
         }
      }

      // Check if both sets contain the same data
      if(std::equal(simSet.begin(), simSet.end(), remSet.begin()))
      {
         return true;
      } else
      {
         return false;
      }
   }

}
