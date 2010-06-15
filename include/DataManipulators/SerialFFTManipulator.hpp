/** \file SerialFFTManipulator.hpp
 *  \brief Serial data manipulator for FFT transforms
 */

#ifndef SERIALFFTMANIPULATOR_HPP
#define SERIALFFTMANIPULATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FFTFlatScalar.hpp"
#include "GeneralScalars/FFTOrderScalar.hpp"
#include "DataManipulators/SerialManipulatorBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Manipulation routines for RTP data
    */
   class SerialFFTManipulator: public SerialManipulatorBase<FFTFlatScalar, FFTOrderScalar>
   {
      public:
         /**
          * @brief Constructor
          *
          * \param pTrunc Pointer to the truncation
          */
         SerialFFTManipulator(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~SerialFFTManipulator() {};

         /**
          * @brief Send data to be converted
          *
          * \param rData Data to send
          */
         void receive(FFTFlatScalar &rData);

         /**
          * @brief Send data to be converted
          *
          * \param rData Data to send
          */
         void receive(FFTOrderScalar &rData);
         
      protected:

      private:
   };
}

#endif // SERIALFFTMANIPULATOR_HPP
