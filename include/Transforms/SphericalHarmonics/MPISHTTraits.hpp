/** \file MPISHTTraits.hpp
 *  \brief Traits for the MPI spherical harmonic transforms
 */

#ifndef MPISHTTRAITS_HPP
#define MPISHTTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GeneralScalars/FDSHOrderScalar.hpp"
#include "GeneralScalars/FFTFlatScalar.hpp"
#include "GeneralScalars/FFTOrderScalar.hpp"
#include "DataManipulators/MPIFFTManipulator.hpp"
#include "Transforms/FFT/FFTransform.hpp"
#include "Transforms/Polynomial/LegendreTransform.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the MPI spherical harmonic transforms
    */
   class MPISHTTraits
   {
      public:
         /// Typedef for the Spherical harmonics scalar data type
         typedef FDSHOrderScalar SHScalarType;

         /// Typedef for the FFT forward data type scalar data type
         typedef FFTFlatScalar FFTForwardType;

         /// Typedef for the FFT backward data type
         typedef FFTOrderScalar FFTBackwardType;

         /// Typedef for the DataManipulator data type
         typedef MPIFFTManipulator DataManipulatorType;

         /// Typedef for the Legendre Transform data type
         typedef LegendreTransform LegendreTransformType;

         /// Typedef for the FFT Transform data type
         typedef FFTransform FFTransformType;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         MPISHTTraits() {};

         /**
          * @brief Destructor
          */
         virtual ~MPISHTTraits() {};
   };

}

#endif // MPISHTTRAITS_HPP
