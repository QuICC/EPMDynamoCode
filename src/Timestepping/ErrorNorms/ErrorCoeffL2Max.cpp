/** \file ErrorCoeffL2Max.cpp
 *  \brief Implementation of the relative L2 max error norm per radial coefficient computation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ErrorNorms/ErrorCoeffL2Max.hpp"

// Project includes
//

namespace EPMDynamo {

   EPMFloat ErrorCoeffL2Max::computeNorm(const ErrorCoeffL2Max::ScalarType& rVar, const ErrorCoeffL2Max::ScalarType& rRef)
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

   ErrorCoeffL2Max::ErrorCoeffL2Max()
   {
   }

}
