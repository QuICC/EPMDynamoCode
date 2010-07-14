/** \file MPIFDSHManipulator.hpp
 *  \brief Implementation of the MPI FDSH data manipulator
 */

#ifndef MPIFDSHMANIPULATOR_HPP
#define MPIFDSHMANIPULATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FDSHOrderScalar.hpp"
#include "GeneralScalars/FDSHDegreeScalar.hpp"
#include "DataManipulators/MPIManipulatorBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Manipulation routines for FDSH data
    *
    * \tparam TForward Type of the forward data
    */
   template <typename TForward> class MPIFDSHManipulator: public MPIManipulatorBase<TForward, FDSHDegreeScalar>
   {
      public:
         /**
          * @brief Constructor
          *
          * \param pTrunc Pointer to the truncation
          */
         MPIFDSHManipulator(const SmartTruncation pTrunc);

         /**
          * @brief Empty destructor
          */
         virtual ~MPIFDSHManipulator() {};

      protected:

      private:
         /**
          * @brief Initialise the MPI datatypes
          */
         void initTypes();

         /**
          * @brief Build an TForward MPI Datatype
          *
          * \param data Input forward data type
          * \param type Created MPI data type
          * \param coreID ID of the cpu/core
          */
         void buildFType(TForward &data, MPI_Datatype &type, const int coreID);

         /**
          * @brief Build an FDSHDegreeScalar MPI Datatype
          *
          * \param data Input data type
          * \param type Created MPI data type
          * \param coreID ID of the cpu/core
          */
         void buildBType(FDSHDegreeScalar &data, MPI_Datatype &type, const int coreID);
   };

   template <typename TForward> MPIFDSHManipulator<TForward>::MPIFDSHManipulator(const SmartTruncation pTrunc)
      : MPIManipulatorBase<TForward, FDSHDegreeScalar>(pTrunc, 2, 2)
   {
      // initialise the MPI datatypes
      this->initTypes();

      // initialise the Buffers sizes
      this->initBufferSizes();
   }

   template <typename TForward> void MPIFDSHManipulator<TForward>::initTypes()
   {
      int nCore = this->mpTrunc->para().nCore();

      // Loop over all cpu
      for(int i = 0; i < nCore; ++i)
      {
         // Create the Forward transform MPI datatypes
         std::map<TForward *, MPI_Datatype>  tFMap;
         for(int j = 0; j < this->mNFTmp; ++j)
         {
            MPI_Datatype   type;
            this->buildFType(this->mFTmp.at(j), type, i);
            tFMap.insert(std::make_pair(&this->mFTmp.at(j), type));
         }
         this->mFTypes.push_back(tFMap);

         // Create the backward transform MPI datatypes
         std::map<FDSHDegreeScalar *, MPI_Datatype>  tBMap;
         for(int j = 0; j < this->mNBTmp; ++j)
         {
            MPI_Datatype   type;
            this->buildBType(this->mBTmp.at(j), type, i);
            tBMap.insert(std::make_pair(&this->mBTmp.at(j), type));
         }
         this->mBTypes.push_back(tBMap);
      }
   }


   template <typename TForward> void MPIFDSHManipulator<TForward>::buildFType(TForward &data, MPI_Datatype &type, const int coreID)
   {
      // Create unique ID out of indexes
      const double nFactor = 1.0e8;
      const double mFactor = 1.0e4;
      const double lFactor = 1.0;

      // Create  map to sort out the required indexes
      std::map<double, ArrayI>  idxSorter;
      // List of local indexes
      std::set<double>  sorterL;
      // List of remote indexes
      std::set<double>  sorterM;
      // List of index in common
      std::set<double>  sorted;
      // Storage for the index
      ArrayI   pos(3);
      
      int n_, m_, l_;
      double key;

      // Create the list of local indexes
      int r0;
      int nR;
      ArrayI   localMs = this->mpTrunc->local()->fdsh()->mArray();
      ArrayI   nL = this->mpTrunc->local()->fdsh()->nLArray(this->mpTrunc->sim()->hoz()->nL());
      for(int m=0; m < localMs.size(); ++m)
      {
         m_ = localMs(m);
         r0 = this->mpTrunc->local()->fdsh()->r0(m);
         nR = this->mpTrunc->local()->fdsh()->nR(m);
         for(int n=0; n < nR; ++n)
         {
            n_ = n + r0;
            for(int l=0; l < nL(m); ++l)
            {
               l_ = m_+l;
               pos(0) = n;
               pos(1) = m;
               pos(2) = l;
               key = n_*nFactor + l_*lFactor + m_*mFactor;
               idxSorter.insert(std::make_pair(key, pos));
               sorterL.insert(key);
            }
         }
      }

      // Create the list of remote indexes
      nR = this->mpTrunc->sim()->rad()->nR();
      ArrayI   remoteLs = this->mpTrunc->remote(coreID)->spec()->lArray();
      for(int l=0; l < this->mpTrunc->remote(coreID)->spec()->nL(); ++l)
      {
         l_ = remoteLs(l);
         for(int m=0; m < this->mpTrunc->remote(coreID)->spec()->nM(l); ++m)
         {
            m_ = this->mpTrunc->remote(coreID)->spec()->m(m,l);
            for(int n=0; n < nR; ++n)
            {
               n_ = n;
               key = n_*nFactor + l_*lFactor + m_*mFactor;
               sorterM.insert(key);
            }
         }
      }

      // Create the list of common indexes
      std::set_intersection(sorterL.begin(), sorterL.end(), sorterM.begin(), sorterM.end(), std::inserter(sorted, sorted.begin()));

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
      int n, l, m;
      int mp = this->mpTrunc->sim()->hoz()->mp();
      int tot = 0;
      for(it = sorted.begin(); it != sorted.end(); ++it)
      {
         n = idxSorter[(*it)](0);
         m = idxSorter[(*it)](1);
         l = idxSorter[(*it)](2);
         m_ = m*mp;
         MPI_Get_address(&data.rCoeff(l,m,n), &element);
         displ[tot] = element - base;
         blocks[tot] = 1;
         tot++;
      }

      // Create MPI datatype
      MPI_Type_create_hindexed(nElements, blocks, displ, MPI_DOUBLE_COMPLEX, &type);
      // Commit MPI datatype
      MPI_Type_commit(&type);
   }

   template <typename TForward> void MPIFDSHManipulator<TForward>::buildBType(FDSHDegreeScalar &data, MPI_Datatype &type, const int coreID)
   {
      // Create unique ID out of indexes
      const double nFactor = 1.0e8;
      const double mFactor = 1.0e4;
      const double lFactor = 1.0;

      // Create  map to sort out the required indexes
      std::map<double, ArrayI>  idxSorter;
      std::set<double>  sorterL;
      std::set<double>  sorterM;
      std::set<double>  sorted;
      ArrayI   pos(3);
      
      // Compute the total number of elements
      int n_, m_, l_;
      double key;

      // Create the list of local indexes
      int nR = this->mpTrunc->sim()->rad()->nR();
      ArrayI   remoteLs = this->mpTrunc->local()->spec()->lArray();
      for(int l=0; l < this->mpTrunc->local()->spec()->nL(); ++l)
      {
         l_ = remoteLs(l);
         for(int m=0; m < this->mpTrunc->local()->spec()->nM(l); ++m)
         {
            m_ = this->mpTrunc->local()->spec()->m(m,l);
            for(int n=0; n < nR; ++n)
            {
               n_ = n;
               pos(0) = n;
               pos(1) = m;
               pos(2) = l;
               key = n_*nFactor + l_*lFactor + m_*mFactor;
               idxSorter.insert(std::pair<double, ArrayI>(key, pos));
               sorterL.insert(key);
            }
         }
      }

      // Create the list of remote indexes
      int r0;
      ArrayI   localMs = this->mpTrunc->remote(coreID)->fdsh()->mArray();
      ArrayI   nL = this->mpTrunc->remote(coreID)->fdsh()->nLArray(this->mpTrunc->sim()->hoz()->nL());
      for(int m=0; m < localMs.size(); ++m)
      {
         m_ = localMs(m);
         r0 = this->mpTrunc->remote(coreID)->fdsh()->r0(m);
         nR = this->mpTrunc->remote(coreID)->fdsh()->nR(m);
         for(int n=0; n < nR; ++n)
         {
            n_ = n + r0;
            for(int l=0; l < nL(m); ++l)
            {
               l_ = m_+l;
               key = n_*nFactor + l_*lFactor + m_*mFactor;
               sorterM.insert(key);
            }
         }
      }

      // Create the list of common indexes
      std::set_intersection(sorterL.begin(), sorterL.end(), sorterM.begin(), sorterM.end(), std::inserter(sorted, sorted.begin()));

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
      int n, l, m;
      int tot = 0;
      for(it = sorted.begin(); it != sorted.end(); ++it)
      {
         n = idxSorter[(*it)](0);
         m = idxSorter[(*it)](1);
         l = idxSorter[(*it)](2);
         MPI_Get_address(&data.rLShell(l)(n, m), &element);
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

#endif // MPIFDSHMANIPULATOR_HPP
