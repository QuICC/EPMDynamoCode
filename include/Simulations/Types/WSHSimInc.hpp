/** \file WSHSimInc.hpp
 *  \brief Simple file to help with inclusion of the right headers
 */

#ifndef WSHSIMINC_HPP
#define WSHSIMINC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GeneralScalars/SpectralSHScalar.hpp"
#include "Polynomials/WorlandPolynomial.hpp"
#include "Operators/DenseOperator.hpp"
#include "Operators/DenseLUOperator.hpp"

#include "Timestepping/IterativeScheme.hpp"
#include "Timestepping/PredictorCorrector/PCTimestepControl.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaMethodImplementation.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluenceMethodImplementation.hpp"

namespace EPMDynamo {
}

#endif // WSHSIMINC_HPP
