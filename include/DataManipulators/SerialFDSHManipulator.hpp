/** \file SerialFDSHManipulator.hpp
 *  \brief Serial data manipulator for FDSH transforms
 */

#ifndef SERIALFDSHMANIPULATOR_HPP
#define SERIALFDSHMANIPULATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FDSHOrderScalar.hpp"
#include "GeneralScalars/FDSHDegreeScalar.hpp"
#include "DataManipulators/SerialManipulatorBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Manipulation routines for FDSH data
    *
    * \tparam TForward Type of data for the forward transform
    */
   template <typename TForward> class SerialFDSHManipulator: public SerialManipulatorBase<TForward, FDSHDegreeScalar>
   {
      public:
         /**
          * @brief Constructor
          *
          * \param pTrunc Pointer to the truncation
          */
         SerialFDSHManipulator(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~SerialFDSHManipulator() {};

         /**
          * @brief Send data
          *
          * \param rData Data to send
          *
          * \epmBug Removed minL(), might be bad idea
          * \epmBug Removed minL(), might be bad idea
          * \epmBug Doesn't work with splitted SH only
          */
         void receive(TForward &rData);

         /**
          * @brief Send data
          *
          * \param rData Data to send
          *
          * \epmBug Removed minL(), might be bad idea
          * \epmBug Doesn't work with splitted SH only
          */
         void receive(FDSHDegreeScalar &rData);
         
      protected:

      private:
   };

   template <typename TForward> SerialFDSHManipulator<TForward>::SerialFDSHManipulator(const SmartTruncation pTrunc)
      : SerialManipulatorBase<TForward, FDSHDegreeScalar>(pTrunc, 3, 3)
   {
   }

   template <typename TForward> void SerialFDSHManipulator<TForward>::receive(TForward &rData)
   {
      // Get data from Send queue
      FDSHDegreeScalar &inLData = this->receiveBData();

      int mp = inLData.mp();
      int m_;
      int nL = inLData.nL();
      int nR = inLData.nR();

      for(int l = 0; l < nL; ++l)
      {
         for(int m = 0; m < inLData.nM(l); ++m)
         {
            m_ = m*mp;
            for(int n=0; n < nR; ++n)
            {
               rData.rCoeff(l-m_,m,n) = inLData.lshell(l)(n,m);
            }
         }
      }

      // Free Temporary storage
      this->freeTmp(inLData);
   }

   template <typename TForward> void SerialFDSHManipulator<TForward>::receive(FDSHDegreeScalar &rData)
   {
      // Get data from Send queue
      TForward  &inMData = this->receiveFData();

      int mp = rData.mp();
      int m_;
      int nL = rData.nL();
      int nR = rData.nR();

      for(int l = 0; l < nL; ++l)
      {
         for(int m = 0; m < rData.nM(l); ++m)
         {
            m_ = m*mp;
            for(int n=0; n < nR; ++n)
            {
               rData.rLShell(l)(n,m) = inMData.coeff(l-m_, m, n);
            }
         }
      }

      // Free Temporary storage
      this->freeTmp(inMData);
   }
}

#endif // SERIALFDSHMANIPULATOR_HPP
