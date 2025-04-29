/** \file SHSplitting.hpp
 *  \brief Implementation of a spherical harmonics load splitting algorithm
 */

#ifndef SHSPLITTING_HPP
#define SHSPLITTING_HPP

// System includes
//
#include <vector>

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
    * \brief Implementation of a spherical harmonics load splitting algorithm
    *
    * \brief Improve load inbalance when doing M splitting
    */
   class SHSplitting: public LoadSplitterBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pSimTrunc Simulation wide truncation information
          * @param nCore Number of cores used
          * @param id ID of the CPU/Core
          */
         SHSplitting(SmartSimTrunc pSimTrunc, const int nCore, const int id);

         /**
          * @brief Destructor
          */
         virtual ~SHSplitting() {};

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
          * @param rR0 Array of start indexes for radial points
          * @param rNr Array of number of radial points
          * @param rMs Array of harmonic orders
          * @param id CPU/Core id
          */
         virtual void splitFDSH(ArrayI &rR0, ArrayI &rNr, ArrayI &rMs, const int id) const;

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
          * @brief Check the splitting limits of algorithm
          */
         void checkLimits() const;

         /**
          * @brief Split the harmonic orders load
          *
          * @param rMs Array of harmonic orders
          * @param id CPU/Core id
          */
         void splitM(ArrayI &rMs, const int id) const;
   };

}

#endif // SHSPLITTING_HPP
