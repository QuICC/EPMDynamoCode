/** \file ErrorModeL2Summed.cpp
 *  \brief Implementation of the relative L2 max error norm per harmonic mode computation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ErrorNorms/ErrorModeL2Summed.hpp"

// Project includes
//
#include "Timestepping/TimestepConfig.hpp"

namespace EPMDynamo {

   EPMFloat ErrorModeL2Summed::computeNorm(const ErrorModeL2Summed::ScalarType& rVar, const ErrorModeL2Summed::ScalarType& rRef)
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
            // Compute error norm
            tmp = rVar.lshell(l).col(m).real().dot(rVar.lshell(l).col(m).real()) + rVar.lshell(l).col(m).imag().dot(rVar.lshell(l).col(m).imag());
            tmp = std::sqrt(tmp);

            // Compute vector norm
            refVal = rRef.lshell(l).col(m).real().dot(rRef.lshell(l).col(m).real()) + rRef.lshell(l).col(m).imag().dot(rRef.lshell(l).col(m).imag());
            refVal = std::sqrt(refVal);

            // Compute mixed absolute-relative error (|err|/(|v|+eta))
            tmp /= refVal + TimestepConfig::TIMESTEP_RELERROR_SCALING;

            // Get summed norm
            norm += tmp;
         }
      }

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      return norm;
   }

   ErrorModeL2Summed::ErrorModeL2Summed()
   {
   }

}
