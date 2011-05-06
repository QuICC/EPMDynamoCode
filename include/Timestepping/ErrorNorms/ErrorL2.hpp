/** \file ErrorL2.hpp
 *  \brief Implementation of the L2 error norm computation
 */

#ifndef ERRORL2_HPP
#define ERRORL2_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Implementation of the L2 error norm computation
    *
    */
   class ErrorL2
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /**
          * @brief Compute the L2 error norm
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
         ErrorL2();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorL2() {};
   };

   inline EPMFloat ErrorL2::computeNorm(const ErrorL2::ScalarType& rVar, const ErrorL2::ScalarType& rRef)
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
            tmp += rVar.lshell(l).col(m).real().dot(rVar.lshell(l).col(m).real()) + rVar.lshell(l).col(m).imag().dot(rVar.lshell(l).col(m).imag());

            // Compute vector norm
            refVal += rRef.lshell(l).col(m).real().dot(rRef.lshell(l).col(m).real()) + rRef.lshell(l).col(m).imag().dot(rRef.lshell(l).col(m).imag());
         }
      }
      tmp = std::sqrt(tmp);
      refVal = std::sqrt(refVal);

      // Compute mixed absolute-relative error (|err|/(|v|+eta))
      tmp /= refVal + TimestepConfig::TIMESTEP_RELERROR_SCALING;

      // Get maximum norm
      norm = std::max(norm, tmp);

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      return norm;
   }

   inline EPMFloat ErrorL2::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

}

#endif // ERRORL2_HPP
