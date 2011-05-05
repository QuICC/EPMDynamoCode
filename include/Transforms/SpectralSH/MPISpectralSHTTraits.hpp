/** \file MPISpectralSHTTraits.hpp
 *  \brief Traits for a MPI radial transform 
 */

#ifndef MPISPECTRALSHTTRAITS_HPP
#define MPISPECTRALSHTTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GeneralFields/TorPolField.hpp"
#include "GeneralScalars/FDSHOrderScalar.hpp"
#include "GeneralScalars/FDSHDegreeScalar.hpp"
#include "DataManipulators/MPIFDSHManipulator.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for a MPI radial transform 
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class MPISpectralSHTTraits
   {
      public:
         /// Typedef for the Spectral Toroidal/Poloidal data type
         typedef TorPolField<TSimType> TorPolFieldType;

         /// Typedef for the FDSH forward data type scalar data type
         typedef FDSHOrderScalar FDSHForwardType;

         /// Typedef for the FDSH backward data type
         typedef FDSHDegreeScalar FDSHBackwardType;

         /// Typedef for the DataManipulator data type
         typedef MPIFDSHManipulator<FDSHForwardType> DataManipulatorType;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         MPISpectralSHTTraits() {};

         /**
          * @brief Destructor
          */
         virtual ~MPISpectralSHTTraits() {};
   };

}

#endif // MPISPECTRALSHTTRAITS_HPP
