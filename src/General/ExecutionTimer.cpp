/** \file ExecutionTimer.cpp
 *  \brief Source of the implementation of a simple run time control
 */

// System includes
//
#include <iostream>

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
      std::cout << "***************** Execution time information *****************" << std::endl;

      std::cout << "     Construction time: " << static_cast<int>(std::ceil(this->mMeanTimes(1))) << " / " << static_cast<int>(std::ceil(this->mMaxTimes(1))) << " / " << static_cast<int>(std::ceil(this->mMinTimes(1))) << " seconds" << std::endl;

      std::cout << "     Initialisation time: " << static_cast<int>(std::ceil(this->mMeanTimes(2))) << " / " << static_cast<int>(std::ceil(this->mMaxTimes(2))) << " / " << static_cast<int>(std::ceil(this->mMinTimes(2))) << " seconds" << std::endl;

      std::cout << "     PreRun time: " << static_cast<int>(std::ceil(this->mMeanTimes(3))) << " / " << static_cast<int>(std::ceil(this->mMaxTimes(3))) << " / " << static_cast<int>(std::ceil(this->mMinTimes(3))) << " seconds" << std::endl;

      std::cout << "     Computation time: " << static_cast<int>(std::ceil(this->mMeanTimes(4))) << " / " << static_cast<int>(std::ceil(this->mMaxTimes(4))) << " / " << static_cast<int>(std::ceil(this->mMinTimes(4))) << " seconds" << std::endl;

      std::cout << "     PostRun time: " << static_cast<int>(std::ceil(this->mMeanTimes(5))) << " / " << static_cast<int>(std::ceil(this->mMaxTimes(5))) << " / " << static_cast<int>(std::ceil(this->mMinTimes(5))) << " seconds" << std::endl;

      std::cout << "--------------------------------------------------------------" << std::endl;

      std::cout << "     Total execution time: " << static_cast<int>(std::ceil(this->mMeanTimes(0))) << " / " << static_cast<int>(std::ceil(this->mMaxTimes(0))) << " / " << static_cast<int>(std::ceil(this->mMinTimes(0))) << " seconds" << std::endl;

      std::cout << "**************************************************************" << std::endl;
   }

}
