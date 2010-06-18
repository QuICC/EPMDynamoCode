/** \file LoadSplitterBase.hpp
 *  \brief Base of the implementation of the load splitting algorithms
 */

#ifndef LOADSPLITTERBASE_HPP
#define LOADSPLITTERBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/SimulationTruncation.hpp"
#include "Domain/CoreTruncation.hpp"
#include "Domain/LocalTruncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the load splitting algorithms
    */
   class LoadSplitterBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pSimTrunc Simulation wide truncation information
          * @param nCore Number of cores used
          * @param id ID of the core
          */
         LoadSplitterBase(SmartSimTrunc pSimTrunc, const int nCore, const int id);

         /**
          * @brief Destructor
          */
         virtual ~LoadSplitterBase() {};
         
      protected:
         /**
          * @brief Get Simulation truncation
          */
         SmartSimTrunc  sim() const;

         /**
          * @brief Get Number of cores
          */
         int nCore() const;

         /**
          * @brief Get id of core
          */
         int id() const;

         /**
          * @brief Split load for the RTP part
          *
          * @param rR0 Start index of radial points
          * @param rNr Number of radial points
          * @param rTh0 Array of start indexes of theta points
          * @param rNth Array of number of theta points
          * @param id CPU/Core id
          */
         virtual void splitRTP(int &rR0, int &rNr, ArrayI &rTh0, ArrayI &rNth, const int id) const = 0;

         /**
          * @brief Split load for the FDSH part
          *
          * @param rR0 Array of start indexes for radial points
          * @param rNr Array of number of radial points
          * @param rMArray Array of harmonic orders
          * @param id CPU/Core id
          */
         virtual void splitFDSH(ArrayI &rR0, ArrayI &rNr, ArrayI &rMArray, const int id) const = 0;

         /**
          * @brief Split load for the Spec part
          *
          * @param rLs Array of harmonic degrees
          * @param rMs Vector of arrays of harmonic orders
          * @param id CPU/Core id
          */
         virtual void splitSpec(ArrayI &rLs, std::vector<ArrayI> &rMs, const int id) const = 0;

         /**
          * @brief Compute a simple balanced split of the elements with regard to the local CPU id
          *
          * @param n0 Start index
          * @param nN Number of indexes
          * @param tot Total number of indexes
          * @param id "Splitting" id
          */
         void balancedSplit(int &n0, int &nN, const int tot, const int id) const;

         /**
          * @brief Compute a simple balanced split of the elements with regard to the given parts
          *
          * @param n0 Start index
          * @param nN Number of indexes
          * @param tot Total number of indexes
          * @param parts "ID" of the splitter
          * @param id "Splitting" id
          */
         void balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id) const;

         /**
          * @brief Compute RTP load sharing
          *
          * @param pLocal Local truncation information
          */
         double computeRTPInbalance(SmartLocalTrunc pLocal) const;

         /**
          * @brief Compute FDSH load sharing
          *
          * @param pLocal Local truncation information
          */
         double computeFDSHInbalance(SmartLocalTrunc pLocal) const;

         /**
          * @brief Compute Spec load sharing
          *
          * @param pLocal Local truncation information
          */
         double computeSpecInbalance(SmartLocalTrunc pLocal) const;

         /**
          * @brief Check RTP splitting consistency
          *
          * @param pRemote Remote truncation information
          */
         bool checkRTPConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const;

         /**
          * @brief Check FDSH splitting consistency
          *
          * @param pRemote Remote truncation information
          */
         bool checkFDSHConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const;

         /**
          * @brief Check Spec splitting consistency
          *
          * @param pRemote Remote truncation information
          */
         bool checkSpecConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const;

      private:
         /**
          * @brief Store Simulation truncation
          */
         SmartSimTrunc  mpSimTrunc;

         /**
          * @brief Number of cores
          */
         int mNcore;

         /**
          * @brief Id of core
          */
         int mID;
   };

   inline SmartSimTrunc LoadSplitterBase::sim() const
   {
      return this->mpSimTrunc;
   }

   inline int LoadSplitterBase::nCore() const
   {
      return this->mNcore;
   }

   inline int LoadSplitterBase::id() const
   {
      return this->mID;
   }

}

#endif // LOADSPLITTERBASE_HPP
