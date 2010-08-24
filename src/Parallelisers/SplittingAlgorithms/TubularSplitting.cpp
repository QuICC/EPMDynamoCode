/** \file TubularSplitting.cpp
 *  \brief Source of the implementation of the massive tubular shaped load splitting algorithm
 */

// System includes
//
#include <map>
#include <queue>
#include <set>
#include <map>

// External includes
//

// Class include
//
#include "Parallelisers/SplittingAlgorithms/TubularSplitting.hpp"

// Project includes
//
#include "General/EPMException.hpp"

namespace EPMDynamo {

   TubularSplitting::TubularSplitting(SmartSimTrunc pSimTrunc, const int nCore, const int id)
      : LoadSplitterBase(pSimTrunc, nCore, id), mHasLongM(false), mMSplitting(0), mRSplitting(0), mMaxSplitting(0)
   {
      // Factorise number of CPUs
      this->factoriseCPUs();

      // Check and set available splitting factors
      this->setSplitting();

      // Initialise the attributes required for computing the M splitting
      this->initLoadM();

      // Compute the load splitting for all CPUs
      this->computeLoadSplitM();
   }

   void TubularSplitting::factoriseCPUs()
   {
      int factor = static_cast<int>(std::sqrt(this->nCore()));

      // Compute smaller factor
      while(factor > 0)
      {
         if(this->nCore() % factor == 0)
         {
            this->mCPUFactors.insert(factor);
         }
         --factor;
      }
   }

   void TubularSplitting::setSplitting()
   {
      // Create list of possible splittings
      std::vector<int> usableR;
      std::set<int>::reverse_iterator  rit = this->mCPUFactors.rbegin();
      bool splittable;
      while(rit != this->mCPUFactors.rend())
      {
         // Set R splitting factor
         this->mRSplitting = *rit;

         // Set M splitting factor
         this->mMSplitting = this->nCore()/this->mRSplitting;

         // Check splitting factors
         splittable = this->checkSplitting();

         if(splittable)
         {
            usableR.push_back(this->mRSplitting);
         }

         // Set M splitting factor
         this->mMSplitting = *rit;

         // Set R splitting factor
         this->mRSplitting = this->nCore()/this->mMSplitting;

         // Check splitting factors
         splittable = this->checkSplitting();

         if(splittable)
         {
            usableR.push_back(this->mRSplitting);
         }

         ++rit;
      }

      // Throw exception if no possibility has been found
      if(usableR.size() < 1)
      {
         throw EPMException("TubularSplitting::checkSplitting", "The problem size and number of CPUs are incompatible!");
      }

      // Choose the best splitting option
      this->chooseBestSplitting(usableR);
   }

   bool TubularSplitting::checkSplitting() const
   {
      // We don't want the extrem cases (no splitting in one direction)
      if( (this->mRSplitting == 1 || this->mMSplitting == 1) && this->nCore() > 2)
      {
         return false;
      }

      // Both splitting factors have to be even (always two in communication)
      if(this->mRSplitting % 2 == 1 || this->mMSplitting % 2 == 1)
      {
         return false;
      }

      // Load balancing requires at least 2 harmonic orders per CPU
      if(this->sim()->hoz()->nM() / this->mMSplitting < 2)
      {
         return false;
      }

      // Load balancing requires that there is at least 1 radial point per CPU
      if(this->sim()->rad()->nR() / this->mRSplitting < 1)
      {
         return false;
      }

      // In all other cases accept the splitting
      return true;
   }

   void TubularSplitting::chooseBestSplitting(const std::vector<int> &usableR)
   {
      // Set M splitting factor
      this->mRSplitting = usableR.at(0);

      // Set M splitting factor
      this->mMSplitting = this->nCore()/this->mRSplitting;

      // Set Max splitting factor
      if(this->mMSplitting >= this->mRSplitting)
      {
         this->mMaxSplitting = this->mMSplitting;
         this->mHasLongM = true;
      } else
      {
         this->mMaxSplitting = this->mRSplitting;
         this->mHasLongM = false;
      }
   }

   void TubularSplitting::initLoadM()
   {
      // Initialise the load list and compute total load
      double totalLoad = 0.0;
      int nL = this->sim()->hoz()->nL();
      int nM = this->sim()->hoz()->nM();
      for(int i = 0; i < nM; ++i)
      {
         this->mLoadListM.push_back(nL - i);
         totalLoad += static_cast<double>(nL -i);
      }

      // Create the list of optimal loads per CPU
      double opt = totalLoad / static_cast<double>(this->mMSplitting);
      int addOpt = static_cast<int>(totalLoad) % this->mMSplitting;
      for(int i = 0; i < this->mMSplitting; ++i)
      {
         this->mLoadOptimalM.push(static_cast<int>(opt) + (i < addOpt));
      }

      // Initialise the per CPU load sums
      for(int i = 0; i < this->mMSplitting; ++i)
      {
         this->mLoadSumM.push_back(0);
      }
   }

   void TubularSplitting::updateLoadSumM()
   {
      typedef  std::multimap<int, int>::iterator MapIt;
      MapIt it;
      std::pair<MapIt, MapIt> range;

      // Loop over number of CPUs
      for(int i = 0; i < this->mMSplitting; ++i)
      {
         // Reset previous sum
         this->mLoadSumM.at(i) = 0;

         // Get assigned loads
         range = this->mLoadSplitM.equal_range(i);

         // Loop over the assigned modes
         for(it = range.first; it != range.second; ++it)
         {
            this->mLoadSumM.at(i) += it->second;
         }
      }
   }

   void TubularSplitting::splitPairsM()
   {
      // Compute the number of rows avaiable per cpus
      int rows = static_cast<int>(std::ceil(static_cast<double>(this->sim()->hoz()->nM())/static_cast<double>(this->mMSplitting)));

      // Check for special cases where care is required in the splitting
      bool oddRows = rows % 2;
      bool needCare = (this->sim()->hoz()->nM() % this->mMSplitting) + oddRows;

      int nPairs;
      int pairMaxM;

      // Need careful splitting ?
      if(needCare)
      {
         // Odd or even number of rows have different algorithm
         if(oddRows)
         {
            nPairs = rows/2 - 1;
         } else
         {
            nPairs = rows/2 - 2;
         }
         pairMaxM = 2*nPairs*this->mMSplitting - 1;

      // Perfect split matching, can simply build all the pairs
      } else
      {
         nPairs = rows/2;

         pairMaxM = this->sim()->hoz()->maxM();
      }
      
      // Counter for the nubmer of elements removed
      int k = 0;

      // Loop over the pairs
      for(int j = 0; j < nPairs; ++j)
      {
         // Loop over the number CPUs
         for(int i = 0; i < this->mMSplitting; ++i)
         {
            this->mLoadSplitM.insert(std::make_pair(i, this->mLoadListM.at(k)));
            this->mLoadSplitM.insert(std::make_pair(i, this->mLoadListM.at(pairMaxM - k)));
            k++;
         }
      }

      // Remove the used loads from the list
      this->mLoadListM.erase(this->mLoadListM.begin(), this->mLoadListM.begin() + 2*k);
   }

   void TubularSplitting::splitRemainingRowsM(const bool hasOddRows)
   {
      // Add a additional row to the assigned loads (in all cases there are at least 2 rows left)
      for(int j = 0; j < this->mMSplitting; ++j)
      {
         this->mLoadSplitM.insert(std::make_pair(j, this->mLoadListM.at(j)));
      }

      // Remove just added loads from list
      this->mLoadListM.erase(this->mLoadListM.begin(), this->mLoadListM.begin() + this->mMSplitting);

      // Special treatment is required if the number of rows is even
      if(! hasOddRows)
      {
         // Additional shifted row (the idea is to get a similar outcome than in the odd case)
         //  Get "upper" half of the cpus
         int uCPU = this->mMSplitting/2 + (this->mMSplitting % 2);
         //  Get "lower" half of the cpus
         int lCPU = this->mMSplitting - uCPU;

         // Assign loads to upper part
         for(int j = 0; j < uCPU; ++j)
         {
            this->mLoadSplitM.insert(std::make_pair(j, this->mLoadListM.at(this->mMSplitting -1  - 2*j)));
         }

         // Assign loads to lower part
         for(int j = 0; j < lCPU; ++j)
         {
            this->mLoadSplitM.insert(std::make_pair(uCPU + j, this->mLoadListM.at(this->mMSplitting - 2  - 2*j)));
         }

         // Remove just added loads from list
         this->mLoadListM.erase(this->mLoadListM.begin(), this->mLoadListM.begin() + this->mMSplitting);
      }

      // Update the load sums per CPU
      this->updateLoadSumM();

      // Distribute the remaining loads
      unsigned int leftLoads = this->mLoadListM.size();
      
      // Progression counters
      unsigned int idx = 0;
      int current;
      int curSum;

      // Algorithm is not perfect, so add margin to fit
      int margin = 0;

      // Loop until there is no remaining load
      while(idx < leftLoads)
      {
         // Reset selected CPU
         current = -1;
         curSum = -1;

         // Search for best place to assign load
         for(int j = 0; j < this->mMSplitting; ++j)
         {
            // Check if position is suitable
            if(this->mLoadSumM.at(j) +  this->mLoadListM.at(idx) <= this->mLoadOptimalM.front() + margin && (this->mLoadSumM.at(j) +  this->mLoadListM.at(idx) > curSum))
            {
               // Set current best fit information
               current = j;
               curSum = this->mLoadSumM.at(j) +  this->mLoadListM.at(idx);

               // Reset margin
               margin = 0;
            }
         }

         // If suitable position has be found assign load to it
         if(current != -1)
         {
            // Assign load
            this->mLoadSplitM.insert(std::make_pair(current,this->mLoadListM.at(idx)));
            // Update related sum
            this->mLoadSumM.at(current) += this->mLoadListM.at(idx);
            // increment counter
            idx++;

            // If CPU reach optimal load, remove it from load queue
            if(this->mLoadSumM.at(current) == this->mLoadOptimalM.front())
            {
               this->mLoadOptimalM.pop();
            }
         // If no suitable position has be found, increase fit margin
         } else
         {
            margin++;
         }
      }

   }

   void TubularSplitting::reshuffleLoadSplitM()
   {
   }

   void TubularSplitting::computeLoadSplitM()
   {
      // Split load in possible pairs
      this->splitPairsM();

      // If required split the remaining loads
      if(this->mLoadListM.size() != 0)
      {
         int rows = static_cast<int>(std::ceil(static_cast<double>(this->mLoadListM.size())/static_cast<double>(this->mMSplitting)));
         bool oddRows =  rows % 2;

         // Split up the remaining rows
         this->splitRemainingRowsM(oddRows);

         // Reshuffle to improve memory balancing
         this->reshuffleLoadSplitM();
      }

      // Update Load sums
      this->updateLoadSumM();
   }

   void TubularSplitting::splitRTP(int &rR0, int &rNr, ArrayI &rTh0, ArrayI &rNth, const int id) const
   {
      // Split the data box in R direction
      this->splitR(rR0, rNr, id);

      // Prepare storage for Theta data box splitting
      rTh0.resize(rNr);
      rNth.resize(rNr);

      int off;
      int n;
      // Split total number of grid points (nR * nTh)
      this->balancedSplit(off, n, this->sim()->hoz()->nTh()*rNr, this->mMSplitting, this->splitMID(id));

      // Compute starting point of grid points
      rTh0.setConstant(0);
      for(int i = 0; i < off; ++i)
      {
         rTh0(i % rNr) += 1;
      }

      // Get smallest th0 to use as offset
      off = 0;
      int minTh0 = rTh0(0);
      for(int r = 1; r < rNr; ++r)
      {
         if(minTh0 > rTh0(r))
         {
            off = r;
            minTh0 = rTh0(r);
         }
      }

      // Compute optimal number of grid points
      rNth.setConstant(0);
      for(int i = 0; i < n; ++i)
      {
         rNth((i + off) % rNr) += 1;
      }
   }

   void TubularSplitting::splitFDSH(ArrayI &rR0, ArrayI &rNr, ArrayI &rMs, const int id) const
   {
      // Split data box in M direction
      this->splitM(rMs, id);

      // Resize the radial splitting informations
      rR0.resize(rMs.size());
      rNr.resize(rMs.size());

      // Split the data box in R direction
      int n;
      int off;
      this->splitR(off, n, id);

      // Set radial sizes and offsets
      rR0.setConstant(off);
      rNr.setConstant(n);
   }

   void TubularSplitting::splitSpec(ArrayI &rLs, std::vector<ArrayI> &rMs, const int id) const
   {
      // clear the vector of array
      rMs.clear();

      // Split data box in M direction
      ArrayI   tmpMs;
      this->splitM(tmpMs, id);

      // Create the queue of all harmonic modes
      std::queue< std::pair<int,int> >  modeQueue;
      int l_;

      // Two degrees get added per iteration
      for(int l = 0; l < this->sim()->hoz()->nL()/2; ++l)
      {
         // Get the 'up-down' degree while l is the 'down-up'
         l_ = this->sim()->hoz()->maxL() - l;

         // Add the 'down-up' mode to queue
         for(int i = 0; i < tmpMs.size(); ++i)
         {
            // Check if harmoni degree and order are compatible
            if(tmpMs(i) <= l)
            {
               modeQueue.push(std::make_pair(l, tmpMs(i)));
            }
         }

         // Add the 'up-down' mode to queue
         for(int i = 0; i < tmpMs.size(); ++i)
         {
            // Check if harmoni degree and order are compatible
            if(tmpMs(i) <= l_)
            {
               modeQueue.push(std::make_pair(l_, tmpMs(i)));
            }
         }
      }

      // Handle the forgotten mode in case of odd total 
      if(this->sim()->hoz()->nL()%2 == 1)
      {
         l_ = this->sim()->hoz()->nL()/2;
         for(int i = 0; i < tmpMs.size(); ++i)
         {
            if(tmpMs(i) <= l_)
            {
               modeQueue.push(std::make_pair(l_, tmpMs(i)));
            }
         }
      }

      int off;
      int nModes;

      // Compute balanced splitting of total number of modes over R factorisation
      this->balancedSplit(off, nModes, modeQueue.size(), this->mRSplitting, this->splitRID(id));

      // Remove offsetted modes from queue
      for(int i = 0; i < off; ++i)
      {
         modeQueue.pop();
      }

      // Prepare sorters for the acquired modes
      std::set<int> sorterL;
      std::multimap<int, int> sorterM;
      std::pair<int, int>  front;

      // Fill in the sorters with the needed modes
      for(int i = 0; i < nModes; ++i)
      {
         front = modeQueue.front();
         modeQueue.pop();
         sorterL.insert(front.first);
         sorterM.insert(front);
      }

      // Resize array of degrees to right size
      rLs.resize(sorterL.size());

      // Store degree in outout array
      std::set<int>::iterator setIt;
      int i = 0;
      for(setIt = sorterL.begin(); setIt != sorterL.end(); ++setIt, ++i)
      {
         rLs(i) = (*setIt);
      }

      // Create vector of arrays of harmonic orders
      std::multimap<int, int>::iterator mmIt;
      std::pair<std::multimap<int, int>::iterator, std::multimap<int, int>::iterator > range;
      for(int j = 0; j < rLs.size(); ++j)
      {
         rMs.push_back(ArrayI(sorterM.count(rLs(j))));
         range = sorterM.equal_range(rLs(j));
         i = 0;
         for(mmIt = range.first; mmIt != range.second; ++mmIt, ++i)
         {
            rMs.at(j)(i) = (*mmIt).second;
         }
      }
   }

   void TubularSplitting::splitR(int &rR0, int &rNr, const int id) const
   {
      this->balancedSplit(rR0, rNr, this->sim()->rad()->nR(), this->mRSplitting, this->splitRID(id));
   }

   void TubularSplitting::splitM(ArrayI &rMs, const int id) const
   {
      typedef  std::multimap<int, int>::const_iterator MapIt;
      MapIt it;
      std::pair<MapIt, MapIt> range;

      // resize the array of harmonic orders
      rMs.resize(this->mLoadSplitM.count(this->splitMID(id)));

      // Get range of harmonic orders for given id
      range = this->mLoadSplitM.equal_range(this->splitMID(id));

      // Put them into set to be sure to get them in a sorted way and convert to orders
      std::set<int>  sorter;
      for(it = range.first; it != range.second; ++it)
      {
         sorter.insert(this->sim()->hoz()->nL() - it->second);
      }

      // Extract the ordered harmonic orders from set and store in output array
      std::set<int>::iterator setIt;
      int i = 0;
      for(setIt = sorter.begin(); setIt != sorter.end(); ++setIt, ++i)
      {
         rMs(i) = *setIt;
      }
   }

   void TubularSplitting::printSpectralInfo() const
   {
      std::cout << "########## IN DEVELOPMENT ###########" << std::endl;
      std::cout << "     Harmonic degree splitting factor: " << this->mRSplitting << std::endl;
      std::cout << "     Harmonic order splitting factor: " << this->mMSplitting << std::endl;
   }

   void TubularSplitting::printSphHarmonicsInfo() const
   {
      std::cout << "########## IN DEVELOPMENT ###########" << std::endl;
      std::cout << "     Radial splitting factor: " << this->mRSplitting << std::endl;
      std::cout << "     Harmonic order splitting factor: " << this->mMSplitting << std::endl;
   }

}
