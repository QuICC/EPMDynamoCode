/** \file ErrorModeL2Max.hpp
 *  \brief Implementation of the relative L2 max error norm per harmonic mode computation
 */

#ifndef ERRORMODEL2MAX_HPP
#define ERRORMODEL2MAX_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the L2 max error norm computation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ErrorModeL2Max
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Compute the Max error norm
          *
          * @param rVar Input variable
          * @param rRef Reference variable for relative error
          */
         static EPMFloat computeNorm(const ScalarType& rVar, const ScalarType& rRef);

         /**
          * @brief Compare new and old values and keep maximum
          *
          * @param newNorm New error norm value
          * @param oldNorm Old error norm value
          */
         static EPMFloat updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm);

      protected:

      private:
         /**
          * @brief Constructor
          */
         ErrorModeL2Max();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorModeL2Max() {}; 
   };

   template <typename TSimType> EPMFloat ErrorModeL2Max<TSimType>::computeNorm(const typename ErrorModeL2Max<TSimType>::ScalarType& rVar, const typename ErrorModeL2Max<TSimType>::ScalarType& rRef)
   {
      // Create temporary storage
      EPMFloat norm = 0.0;
      EPMFloat refVal = 0.0;
      EPMFloat tmp = 0.0;

      // Get truncation information
      const int l0 = rVar.minL();
      int nL = rVar.nL();
      int nM;

      // Loop over harmonic modes
      for(int l = l0; l < nL; ++l)
      {
         nM = rVar.nM(l);
         for(int m=0; m < nM; ++m)
         {
            // Compute L2 norm            
            tmp = rVar.lshell(l).col(m).real().dot(rVar.lshell(l).col(m).real()) + rVar.lshell(l).col(m).imag().dot(rVar.lshell(l).col(m).imag());
            tmp = std::sqrt(tmp);
            refVal = rRef.lshell(l).col(m).real().dot(rRef.lshell(l).col(m).real()) + rRef.lshell(l).col(m).imag().dot(rRef.lshell(l).col(m).imag());
            refVal = std::sqrt(refVal);
            if(refVal > TimestepConfig::TIMESTEP_RELERROR_THRESHOLD)
            {
               tmp /= refVal;
            }
            norm = std::max(norm, tmp);
         }
      }

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      return norm;
   }

   template <typename TSimType> inline EPMFloat ErrorModeL2Max<TSimType>::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

   template <typename TSimType> ErrorModeL2Max<TSimType>::ErrorModeL2Max()
   {
   }

}

#endif // ERRORMODEL2MAX_HPP
