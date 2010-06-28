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

#include "Timestepping/PredictorCorrector/PCScheme.hpp"
#include "Timestepping/PredictorCorrector/PCTimestepControl.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaMethod.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluenceMethod.hpp"

#include "Timestepping/ETD/ETDScheme.hpp"
#include "Timestepping/ETD/ETDTimestepControl.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKMethod.hpp"
#include "Timestepping/ETD/ETD2RK/ETD2RKInfluenceMethod.hpp"

namespace EPMDynamo {
}

#endif // WSHSIMINC_HPP
