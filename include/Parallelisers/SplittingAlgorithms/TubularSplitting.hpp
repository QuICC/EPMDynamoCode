/** \file TubularSplitting.hpp
 *  \brief Implemenation of a massive load splitting in a tubular fashion
 */

#ifndef TUBULARSPLITTING_HPP
#define TUBULARSPLITTING_HPP

// System includes
//
#include <vector>
#include <set>
#include <deque>
#include <queue>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/SimulationTruncation.hpp"
#include "Domain/CoreTruncationBase.hpp"
#include "Domain/CoreTruncation.hpp"
#include "Domain/LocalTruncation.hpp"
#include "Parallelisers/LoadSplitterBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of a massive load splitting in a tubular fashion
    *
    * \epmTodo Review implementation and improve splitting algorithm. The choice of splitting factors can most likely be improved!
    */
   class TubularSplitting: public LoadSplitterBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pSimTrunc Simulation wide truncation information
          * @param nCore Number of cores used
          * @param id ID of the CPU/Core
          */
         TubularSplitting(SmartSimTrunc pSimTrunc, const int nCore, const int id);

         /**
          * @brief Destructor
          */
         virtual ~TubularSplitting() {};

         /**
          * @brief Print splitting information for spectral transform
          */
         void printSpectralInfo() const;

         /**
          * @brief Print splitting information for spherical harmonics transform
          */
         void printSphHarmonicsInfo() const;
         
      protected:

         /**
          * @brief Split load for the RTP part
          *
          * @param rR0 Start index of radial points
          * @param rNr Number of radial points
          * @param rTh0 Array of start indexes of theta points
          * @param rNth Array of number of theta points
          * @param id CPU/Core id
          */
         virtual void splitRTP(int &rR0, int &rNr, ArrayI &rTh0, ArrayI &rNth, const int id) const ;

         /**
          * @brief Split load for the FDSH part
          *
          * @param rLs Array of harmonic degrees
          * @param rM0 Array of start indexes of harmonic orders
          * @param rNm Array of number of harmonic orders
          * @param id CPU/Core id
          */
         virtual void splitFDSH(ArrayI &rLs, ArrayI &rM0, ArrayI &rNm, const int id) const;

         /**
          * @brief Split load for the Spec part
          *
          * @param rLs Array of harmonic degrees
          * @param rMs Vector of arrays of harmonic orders
          * @param id CPU/Core id
          */
         virtual void splitSpec(ArrayI &rLs, std::vector<ArrayI> &rMs, const int id) const;

      private:
         /**
          * @brief Bool for choosing between long or short M
          */
         bool mHasLongM;

         /**
          * @brief Splitting factor for the M dimension
          */
         int mMSplitting;

         /**
          * @brief Splitting factor for the R dimension
          */
         int mRSplitting;

         /**
          * @brief Max Splitting factor (is used for the ordering) 
          */
         int mMaxSplitting;

         /**
          * @brief Factors of the number of CPUs
          */
         std::set<int> mCPUFactors;

         /**
          * @brief Storage for the list of loads
          */
         std::deque<int>   mLoadListM;

         /**
          * @brief Storage for the list of optimal per CPU loads
          */
         std::queue<int>   mLoadOptimalM;

         /**
          * @brief Storage for summed load per CPU
          */
         std::vector<int>   mLoadSumM;

         /**
          * @brief Storage for the splitted load per CPU
          */
         std::multimap<int, int>   mLoadSplitM;

         /**
          * @brief Decompose number of CPUs into factors
          */
         void factoriseCPUs();

         /**
          * @brief Choose the splitting factors
          */
         void setSplitting();

         /**
          * @brief Check the splitting factors compatibility
          */
         bool checkSplitting() const;

         /**
          * @brief Choose best splitting option
          *
          * @param usableR vector of usable ratios
          */
         void chooseBestSplitting(const std::vector<int> &usableR);

         /**
          * @brief Initialise the attributes required for the M splitting
          */
         void initLoadM();

         /**
          * @brief Update the load sum per CPU
          */
         void updateLoadSumM();

         /**
          * @brief Split the possible pairs of M values
          */
         void splitPairsM();

         /**
          * @brief Split the load for the remaining rows after pairing
          *
          * @param hasOddRows Is the number of rows odd?
          */
         void splitRemainingRowsM(const bool hasOddRows);

         /**
          * @brief Precompute the load splitting for M
          */
         void computeLoadSplitM();

         /**
          * @brief Reshuffle splitting to improve memory load balancing
          */
         void reshuffleLoadSplitM();

         /**
          * @brief Split R direction
          *
          * @param rR0 Radial start index
          * @param rNr Number of radial indexes
          * @param id ID of CPU/Core
          */
         void splitR(int &rR0, int &rNr, const int id) const;

         /**
          * @brief Split M direction
          *
          * @param rMs Array of harmonic orders
          * @param id ID of CPU/Core
          */
         void splitM(ArrayI &rMs, const int id) const;

         /**
          * @brief Convert ID to long axis
          *
          * @param id ID of CPU/Core
          */
         int splitMID(const int id) const;

         /**
          * @brief Convert ID to short axis
          *
          * @param id ID of CPU/Core
          */
         int splitRID(const int id) const;
   };

   inline int TubularSplitting::splitMID(const int id) const
   {
      if(this->mHasLongM)
      {
         return id % this->mMaxSplitting;
      } else
      {
         return id / this->mMaxSplitting;
      }
   }

   inline int TubularSplitting::splitRID(const int id) const
   {
      if(this->mHasLongM)
      {
         return id / this->mMaxSplitting;
      } else
      {
         return id % this->mMaxSplitting;
      }
   }

}

#endif // TUBULARSPLITTING_HPP
