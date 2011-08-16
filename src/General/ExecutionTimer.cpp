/** \file ExecutionTimer.cpp
 *  \brief Source of the implementation of a simple run time control
 */

// System includes
//
#include <iostream>
#include <iomanip>

// External includes
//

// Class include
//
#include "General/ExecutionTimer.hpp"

// Project includes
//

namespace EPMDynamo {

   ExecutionTimer::ExecutionTimer(const bool autostart)
      : TimerType(autostart), mMeanTimes(6), mMaxTimes(6), mMinTimes(6)
   {
      // initialise the times
      this->mMeanTimes.setConstant(0.0);
      this->mMaxTimes.setConstant(0.0);
      this->mMinTimes.setConstant(0.0);
   }

   void ExecutionTimer::update(const int id)
   {
      // Increment measured time
      this->mMeanTimes(id) += this->time();

      // Id = 0 simply sums up every thing
      if(id != 0)
      {
         this->mMeanTimes(0) += this->time();
      }
   }

   void ExecutionTimer::analyze()
   {
      // Get the "global" times from MPI code
      #ifdef EPMDYNAMO_MPI
         // Get the max values
         MPI_Allreduce(this->mMeanTimes.data(), this->mMaxTimes.data(), 7, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         // Get the min values
         MPI_Allreduce(this->mMeanTimes.data(), this->mMinTimes.data(), 7, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
         // Get the mean values
         MPI_Allreduce(MPI_IN_PLACE, this->mMeanTimes.data(), 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         int size;
         MPI_Comm_size(MPI_COMM_WORLD, &size);
         this->mMeanTimes /= static_cast<EPMFloat>(size);

      #endif
   }

   void ExecutionTimer::printInfo() const
   {
      std::cout << std::endl;
      std::cout << "************* Execution time information *************" << std::endl;

      std::cout << "     Construction time: " << std::fixed << std::setprecision(1) << this->mMeanTimes(1) << " / " << this->mMaxTimes(1) << " / " << this->mMinTimes(1) << " seconds" << std::endl;

      std::cout << "     Initialisation time: " << this->mMeanTimes(2) << " / " << this->mMaxTimes(2) << " / " << this->mMinTimes(2) << " seconds" << std::endl;

      std::cout << "     PreRun time: " << this->mMeanTimes(3) << " / " << this->mMaxTimes(3) << " / " << this->mMinTimes(3) << " seconds" << std::endl;

      std::cout << "     Computation time: " << this->mMeanTimes(4) << " / " << this->mMaxTimes(4) << " / " << this->mMinTimes(4) << " seconds" << std::endl;

      std::cout << "     PostRun time: " << this->mMeanTimes(5) << " / " << this->mMaxTimes(5) << " / " << this->mMinTimes(5) << " seconds" << std::endl;

      std::cout << "------------------------------------------------------" << std::endl;

      std::cout << "     Total execution time: " << this->mMeanTimes(0) << " / " << this->mMaxTimes(0) << " / " << this->mMinTimes(0) << " seconds" << std::endl;

      std::cout << "******************************************************" << std::endl;
   }

}
