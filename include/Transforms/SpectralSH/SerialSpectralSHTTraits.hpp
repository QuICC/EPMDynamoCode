/** \file SerialSpectralSHTTraits.hpp
 *  \brief Traits for a serial radial transform
 */

#ifndef SERIALSPECTRALSHTTRAITS_HPP
#define SERIALSPECTRALSHTTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GeneralFields/TorPolField.hpp"
#include "GeneralScalars/FDSHOrderScalar.hpp"
#include "GeneralScalars/FDSHDegreeScalar.hpp"
#include "GeneralFields/TorPolField.hpp"
#include "DataManipulators/SerialFDSHManipulator.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for a serial radial transform 
    */
   class SerialSpectralSHTTraits
   {
      public:
         /// Typedef for the Spectral Toroidal/Poloidal data type
         typedef TorPolField TorPolFieldType;

         /// Typedef for the FDSH forward data type scalar data type
         typedef FDSHOrderScalar FDSHForwardType;

         /// Typedef for the FDSH backward data type
         typedef FDSHDegreeScalar FDSHBackwardType;

         /// Typedef for the DataManipulator data type
         typedef SerialFDSHManipulator<FDSHForwardType> DataManipulatorType;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         SerialSpectralSHTTraits() {};

         /**
          * @brief Destructor
          */
         virtual ~SerialSpectralSHTTraits() {};
   };

}

#endif // SERIALSPECTRALSHTTRAITS_HPP
