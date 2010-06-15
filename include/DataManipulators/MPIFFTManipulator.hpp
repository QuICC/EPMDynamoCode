/** \file MPIFFTManipulator.hpp
 *  \brief Implementation of the MPI FFT data manipulator
 */

#ifndef MPIFFTMANIPULATOR_HPP
#define MPIFFTMANIPULATOR_HPP

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
#include "DataManipulators/MPIManipulatorBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Manipulation routines for FFT data
    *
    * \bug Review/Clean code
    */
   class MPIFFTManipulator: public MPIManipulatorBase<FFTFlatScalar, FFTOrderScalar>
   {
      public:
         /**
          * @brief Constructor
          *
          * \param pTrunc Pointer to the truncation
          */
         MPIFFTManipulator(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~MPIFFTManipulator() {};

      protected:

      private:
         /**
          * @brief Initialise the MPI datatypes
          */
         void initTypes();

         /**
          * @brief Build an FFTFlatScalar MPI Datatype
          *
          * \param data Input data
          * \param type Created MPI data type
          * \param coreID ID of the cpu/core
          */
         void buildFType(FFTFlatScalar &data, MPI_Datatype &type, const int coreID);

         /**
          * @brief Build an FFTOrderScalar MPI Datatype
          *
          * \param data Input data
          * \param type Created MPI data type
          * \param coreID ID of the cpu/core
          */
         void buildBType(FFTOrderScalar &data, MPI_Datatype &type, const int coreID);
   };
}

#endif // MPIFFTMANIPULATOR_HPP
