/** \file TransformsInc.hpp
 *  \brief Simple file to help with inclusion of the right headers
 */

#ifndef TRANSFORMSINC_HPP
#define TRANSFORMSINC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Transforms/SpectralSH/TorPolSpectralSHTransform.hpp"
#include "Transforms/SpectralSH/SerialSpectralSHTTraits.hpp"

#ifdef EPMDYNAMO_MPI
   #include "Transforms/SpectralSH/MPISpectralSHTTraits.hpp"
   #include "Transforms/SphericalHarmonics/MPISHTTraits.hpp"
#endif //EPMDYNAMO_MPI
#include "Transforms/SphericalHarmonics/SerialSHTTraits.hpp"

namespace EPMDynamo {
}

#endif // TRANSFORMSINC_HPP
