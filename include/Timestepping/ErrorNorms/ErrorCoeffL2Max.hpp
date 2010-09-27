/** \file ErrorCoeffL2Max.hpp
 *  \brief Implementation of the relative L2 max error norm per radial coefficient computation
 */

#ifndef ERRORCOEFFL2MAX_HPP
#define ERRORCOEFFL2MAX_HPP

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
   template <typename TSimType> class ErrorCoeffL2Max
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
         ErrorCoeffL2Max();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorCoeffL2Max() {}; 
   };

   template <typename TSimType> EPMFloat ErrorCoeffL2Max<TSimType>::computeNorm(const typename ErrorCoeffL2Max<TSimType>::ScalarType& rVar, const typename ErrorCoeffL2Max<TSimType>::ScalarType& rRef)
   {
      // Create temporary storage
      EPMFloat norm = 0.0;
      EPMFloat refVal = 0.0;
      EPMFloat tmp = 0.0;

      // Get truncation information
      const int l0 = rVar.minL();
      int nL = rVar.nL();
      int nN = rVar.nN();
      int nM;

      // Loop over harmonic modes
      for(int l = l0; l < nL; ++l)
      {
         nM = rVar.nM(l);
         for(int m=0; m < nM; ++m)
         {
            // Compute L2 norm
            for(int n=0; n < nN; ++n)
            {
               // Compute error norm
               tmp = rVar.lshell(l)(n,m).real()*rVar.lshell(l)(n,m).real() + rVar.lshell(l)(n,m).imag()*rVar.lshell(l)(n,m).imag();
               tmp = std::sqrt(tmp);

               // Compute vector norm
               refVal = rRef.lshell(l)(n,m).real()*rRef.lshell(l)(n,m).real() + rRef.lshell(l)(n,m).imag()*rRef.lshell(l)(n,m).imag();
               refVal = std::sqrt(refVal);

               // Compute mixed absolute-relative error (|err|/(|v|+eta))
               tmp /= refVal + TimestepConfig::TIMESTEP_RELERROR_SCALING;

               // Get maximum norm
               norm = std::max(norm, tmp);
            }
         }
      }

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      return norm;
   }

   template <typename TSimType> inline EPMFloat ErrorCoeffL2Max<TSimType>::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

   template <typename TSimType> ErrorCoeffL2Max<TSimType>::ErrorCoeffL2Max()
   {
   }

}

#endif // ERRORCOEFFL2MAX_HPP
