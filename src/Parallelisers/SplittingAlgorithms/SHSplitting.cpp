/** \file SHSplitting.cpp
 *  \brief Source of spherical harmonics load splitting algorithm
 */

// System includes
//
#include <utility>
#include <queue>
#include <set>
#include <assert.h>

// External includes
//

// Class include
//
#include "Parallelisers/SplittingAlgorithms/SHSplitting.hpp"

// Project includes
//
#include "General/EPMException.hpp"

namespace EPMDynamo {

   SHSplitting::SHSplitting(SmartSimTrunc pSimTrunc, const int nCore, const int id)
      : LoadSplitterBase(pSimTrunc, nCore, id)
   {
      // Check the limits of the algorithms
      this->checkLimits();
   }

   void SHSplitting::splitRTP(int &rR0, int &rNr, ArrayI &rTh0, ArrayI &rNth, const int id) const
   {
      int off;
      int n;
      const int nTh = this->sim()->hoz()->nTh();

      // Split total number of grid points (nR * nTh)
      this->balancedSplit(off, n, nTh*this->sim()->rad()->nR(), this->nCore(), id);

      // Set the radial truncation information (full truncation)
      rR0 = off/nTh;
      rNr = std::ceil(static_cast<double>(n)/static_cast<double>(nTh));

      // Resize theta truncation array
      rTh0.resize(rNr);
      rTh0.setConstant(0);
      rNth.resize(rNr);

      // Compute first offset
      rTh0(0) = off % nTh;

      rNth.setConstant(0);
      for(int i = 0; i < n; ++i)
      {
         rNth((rTh0(0)+i)/nTh)++; 
      }
   }

   void SHSplitting::splitFDSH(ArrayI &rR0, ArrayI &rNr, ArrayI &rMs, const int id) const
   {
      // Split the harmonic orders over cpu
      this->splitM(rMs, id);

      // Set the radial truncation information (full truncation)
      rR0.resize(rMs.size());
      rR0.setConstant(0);

      rNr.resize(rMs.size());
      rNr.setConstant(this->sim()->rad()->nR());
   }

   void SHSplitting::splitSpec(ArrayI &rLs, std::vector<ArrayI> &rMs, const int id) const
   {
      // temporary storage for the order splitting
      ArrayI tmpMs;

      // Split the harmonic orders over cpu
      this->splitM(tmpMs, id);

      // Set the harmonic degrees truncation information (full truncation)
      rLs = this->sim()->hoz()->lArray().end(this->sim()->hoz()->nL()- tmpMs(0));

      // Set the harmonic orders truncation information
      rMs.clear();
      int nMs;
      for(int l=0; l < rLs.size(); ++l)
      {
         // Count the number of M applicable to l
         nMs = 0;
         for(int i = 0; i < tmpMs.size(); ++i)
         {
            if(tmpMs(i) <= rLs(l))
            {
               nMs++;
            }
         }

         // Create storage for the orders
         rMs.push_back(ArrayI(nMs));

         // Fill in the orders
         for(int j = 0; j < rMs.at(l).size(); ++j)
         {
            rMs.at(l)(j) = tmpMs(j);
         }
      }
   }

   void SHSplitting::splitM(ArrayI &rMs, const int id) const
   {
      // Get the array of all m
      ArrayI m = this->sim()->hoz()->mArray();
      int nM =  this->sim()->hoz()->nM();

      rMs.resize(m.size()/this->nCore() + (id < nM % this->nCore()));

      int j = 0;
      for(int i = 0; i < nM; ++i)
      {
         if(i % this->nCore() == id)
         {
            rMs(j) = m(i);
            j++;
         }
      }
   }

   void SHSplitting::checkLimits() const
   {
      if(this->nCore() > (this->sim()->hoz()->nM()+1)/2)
      {
         throw EPMException("SHSplitting::checkLimits", "The number of CPUs is too big compared to harmonic order truncation!");
      }

      if(this->nCore() > this->sim()->rad()->nR() * this->sim()->hoz()->nTh())
      {
         throw EPMException("SHSplitting::checkLimits", "The number of CPUs is too big compared to the RTP resolution!");
      }
   }

   void SHSplitting::printSpectralInfo() const
   {
      std::cout << "        No load splitting required!" << std::endl;
   }

   void SHSplitting::printSphHarmonicsInfo() const
   {
      std::cout << "########## IN DEVELOPMENT ###########" << std::endl;
   }

}
