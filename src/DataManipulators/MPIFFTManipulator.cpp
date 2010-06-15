/** \file MPIFFTManipulator.cpp
 *  This file contains the definitions of the MPIFFTManipulator class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "DataManipulators/MPIFFTManipulator.hpp"

// Project includes
//

namespace EPMDynamo {

   MPIFFTManipulator::MPIFFTManipulator(const SmartSTrunc pSTrunc)
      : MPIManipulatorBase<FFTFlatScalar, FFTOrderScalar>(pSTrunc, 2, 2)
   {
      // initialise the MPI datatypes
      this->initTypes();

      // initialise the Buffers sizes
      this->initBufferSizes();
   }

   void MPIFFTManipulator::initTypes()
   {
      int nCore = mpSTrunc->para().nCore();
      for(int i = 0; i < nCore; ++i)
      {
         std::map<FFTFlatScalar *, MPI_Datatype>  tFMap;
         for(int j = 0; j < this->mNFTmp; ++j)
         {
            MPI_Datatype   type;
            this->buildFType(this->mFTmp.at(j), type, i);
            tFMap.insert(std::make_pair(&this->mFTmp.at(j), type));
         }
         this->mFTypes.push_back(tFMap);

         std::map<FFTOrderScalar *, MPI_Datatype>  tBMap;
         for(int j = 0; j < mNBTmp; ++j)
         {
            MPI_Datatype   type;
            this->buildBType(this->mBTmp.at(j), type, i);
            tBMap.insert(std::make_pair(&this->mBTmp.at(j), type));
         }
         this->mBTypes.push_back(tBMap);
      }
   }

   void MPIFFTManipulator::buildFType(FFTFlatScalar &data, MPI_Datatype &type, const int coreID)
   {
      // Create unique ID out of indexes
      const double nFactor = 1.0e8;
      const double mFactor = 1.0e4;
      const double thFactor = 1.0;

      // Create  map to sort out the required indexes
      std::map<double, ArrayI>  idxSorter;
      // List of local indexes
      std::set<double>  sorterM;
      // List of remote indexes
      std::set<double>  sorterF;
      // List of index in common
      std::set<double>  sorted;
      // Storage for the index
      ArrayI   pos(3);
      
      int n_, th_, m_;
      double key;

      // Create the list of local indexes
      int r0 = this->mpSTrunc->local()->rtp()->r0();
      int nR = this->mpSTrunc->local()->rtp()->nR();
      ArrayI localNth = this->mpSTrunc->local()->rtp()->nThArray();
      ArrayI localTh0 = this->mpSTrunc->local()->rtp()->th0Array();
      int nM = this->mpSTrunc->sim()->hoz()->nM();

      // Create the list of local indexes
      for(int n=0; n < nR; ++n)
      {
         n_ = n + r0;
         for(int th=0; th < localNth(n); ++th)
         {
            th_ = localTh0(n) + th;
            for(int m=0; m < nM; ++m)
            {
               pos(0) = n;
               pos(1) = th;
               pos(2) = m;
               key = n_*nFactor + th_*thFactor + m*mFactor;
               idxSorter.insert(std::make_pair(key, pos));
               sorterF.insert(key);
            }
         }
      }

      // Create the list of remote indexes
      int nTh = this->mpSTrunc->sim()->hoz()->nTh();
      ArrayI   remoteMs = this->mpSTrunc->remote(coreID)->fdsh()->mArray();
      for(int m=0; m < this->mpSTrunc->remote(coreID)->fdsh()->nM(); ++m)
      {
         m_ = remoteMs(m);
         for(int n=0; n < this->mpSTrunc->remote(coreID)->fdsh()->nR(m); ++n)
         {
            n_ = mpSTrunc->remote(coreID)->fdsh()->r0(m) + n;
            for(int th=0; th < nTh; ++th)
            {
               key = n_*nFactor + th*thFactor + m_*mFactor;
               sorterM.insert(key);
            }
         }
      }

      // Create the list of common indexes
      std::set_intersection(sorterM.begin(), sorterM.end(), sorterF.begin(), sorterF.end(), std::inserter(sorted, sorted.begin()));

      // Get the number of elements
      int nElements = sorted.size();
      MPI_Aint    displ[nElements];
      int         blocks[nElements];

      MPI_Aint    base;
      MPI_Aint    element;
      //MPI_Get_address(&data, &base);
      base = 0;

      std::set<double>::iterator it;

      // Create MPI displacement list
      int n, th, m;
      int tot = 0;
      for(it = sorted.begin(); it != sorted.end(); ++it)
      {
         n = idxSorter[(*it)](0);
         th = idxSorter[(*it)](1);
         m = idxSorter[(*it)](2);
         if(n == 0)
         {
            nTh = 0;
         } else
         {
            nTh = localNth.segment(0,n).sum();
         }
         MPI_Get_address(&data.rFlat()(m, nTh + th), &element);
         displ[tot] = element - base;
         blocks[tot] = 1;
         tot++;
      }

      // Create MPI datatype
      MPI_Type_create_hindexed(nElements, blocks, displ, MPI_DOUBLE_COMPLEX, &type);
      // Commit MPI datatype
      MPI_Type_commit(&type);
   }

   void MPIFFTManipulator::buildBType(FFTOrderScalar &data, MPI_Datatype &type, const int coreID)
   {
      // Create unique ID out of indexes
      const double nFactor = 1.0e8;
      const double mFactor = 1.0e4;
      const double thFactor = 1.0;

      // Create  map to sort out the required indexes
      std::map<double, ArrayI>  idxSorter;
      std::set<double>  sorterM;
      std::set<double>  sorterF;
      std::set<double>  sorted;
      ArrayI   pos(3);
      
      int n_, th_, m_;
      double key;

      // Create the list of local indexes
      int r0;
      int nR;
      int nTh = mpSTrunc->sim()->hoz()->nTh();
      int nM = mpSTrunc->local()->fdsh()->nM();

      // Create the list of local indexes
      ArrayI   localMs = mpSTrunc->local()->fdsh()->mArray();
      for(int m=0; m < nM; ++m)
      {
         r0 = mpSTrunc->local()->fdsh()->r0(m);
         nR = mpSTrunc->local()->fdsh()->nR(m);
         m_ = localMs(m);
         for(int n=0; n < nR; ++n)
         {
            n_ = r0 + n;
            for(int th=0; th < nTh; ++th)
            {
               pos(0) = n;
               pos(1) = th;
               pos(2) = m;
               key = n_*nFactor + th*thFactor + m_*mFactor;
               idxSorter.insert(std::make_pair(key, pos));
               sorterM.insert(key);
            }
         }
      }

      // Create the list of remote indexes
      r0 = mpSTrunc->remote(coreID)->rtp()->r0();
      nR = mpSTrunc->remote(coreID)->rtp()->nR();
      ArrayI remoteNth = mpSTrunc->remote(coreID)->rtp()->nThArray();
      ArrayI remoteTh0 = mpSTrunc->remote(coreID)->rtp()->th0Array();
      nM = mpSTrunc->sim()->hoz()->nM();
      for(int n=0; n < nR; ++n)
      {
         n_ = n + r0;
         for(int th=0; th < remoteNth(n); ++th)
         {
            th_ = remoteTh0(n) + th;
            for(int m=0; m < nM; ++m)
            {
               key = n_*nFactor + th_*thFactor + m*mFactor;
               sorterF.insert(key);
            }
         }
      }

      // Create the list of common indexes
      std::set_intersection(sorterM.begin(), sorterM.end(), sorterF.begin(), sorterF.end(), std::inserter(sorted, sorted.begin()));

      // Get the number of elements
      int nElements = sorted.size();
      MPI_Aint    displ[nElements];
      int         blocks[nElements];

      MPI_Aint    base;
      MPI_Aint    element;
      //MPI_Get_address(&data, &base);
      base = 0;

      // Create MPI displacement list
      std::set<double>::iterator it;
      int n, th, m;
      int tot = 0;
      for(it = sorted.begin(); it != sorted.end(); ++it)
      {
         n = idxSorter[(*it)](0);
         th = idxSorter[(*it)](1);
         m = idxSorter[(*it)](2);
         MPI_Get_address(&data.rMShell(m)(th, n), &element);
         displ[tot] = element - base;
         blocks[tot] = 1;
         tot++;
      }

      // Create MPI datatype
      MPI_Type_create_hindexed(nElements, blocks, displ, MPI_DOUBLE_COMPLEX, &type);
      // Commit MPI datatype
      MPI_Type_commit(&type);
   }

}
