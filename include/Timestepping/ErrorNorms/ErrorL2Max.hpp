/** \file ErrorL2Max.hpp
 *  \brief Implementation of the L2 max error norm computation
 */

#ifndef ERRORL2MAX_HPP
#define ERRORL2MAX_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the L2 max error norm computation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ErrorL2Max
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Compute the Max error norm
          *
          * @param rVar Input variable
          */
         static DynamoFloat computeNorm(const ScalarType& rVar);

         /**
          * @brief Compare new and old values and keep maximum
          *
          * @param newNorm New error norm value
          * @param oldNorm Old error norm value
          */
         static DynamoFloat updateNorm(const DynamoFloat newNorm, const DynamoFloat oldNorm);

      protected:

      private:
         /**
          * @brief Constructor
          */
         ErrorL2Max();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorL2Max() {}; 
   };

   template <typename TSimType> DynamoFloat ErrorL2Max<TSimType>::computeNorm(const typename ErrorL2Max<TSimType>::ScalarType& rVar)
   {
      // Create temporary storage
      DynamoFloat norm = 0.0;
      DynamoFloat tmp = 0.0;

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
            norm = std::max(norm, std::sqrt(tmp));
         }
      }

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      return norm;
   }

   template <typename TSimType> inline DynamoFloat ErrorL2Max<TSimType>::updateNorm(const DynamoFloat newNorm, const DynamoFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

   template <typename TSimType> ErrorL2Max<TSimType>::ErrorL2Max()
   {
   }

}

#endif // ERRORL2MAX_HPP
