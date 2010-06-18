/** \file PolyOpsTypes.hpp
 */

#ifndef POLYOPSTYPES_HPP
#define POLYOPSTYPES_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/PolyOps/POperator.hpp"
#include "Polynomials/PolyOps/LPOperator.hpp"
#include "Polynomials/PolyOps/RPOperator.hpp"
#include "Polynomials/PolyOps/LRPOperator.hpp"
#include "Polynomials/PolyOps/CPOperator.hpp"
#include "Polynomials/PolyOps/CLPOperator.hpp"
#include "Polynomials/PolyOps/CRPOperator.hpp"
#include "Polynomials/PolyOps/CLRPOperator.hpp"

namespace EPMDynamo {

   /// Typedef for a smart POperator
   typedef EPMSHARED_PTR<PolynomialOperator<POperator> >  SmartPOp;

   /// Typedef for a smart LPOperator
   typedef EPMSHARED_PTR<PolynomialOperator<LPOperator> >  SmartLPOp;

   /// Typedef for a smart RPOperator
   typedef EPMSHARED_PTR<PolynomialOperator<RPOperator> >  SmartRPOp;

   /// Typedef for a smart LRPOperator
   typedef EPMSHARED_PTR<PolynomialOperator<LRPOperator> >  SmartLRPOp;

   /// Typedef for a smart CPOperator with real constant
   typedef EPMSHARED_PTR<PolynomialOperator<CPOperator<EPMFloat> > >  SmartDPOp;

   /// Typedef for a smart CLPOperator with real constant
   typedef EPMSHARED_PTR<PolynomialOperator<CLPOperator<EPMFloat> > >  SmartDLPOp;

   /// Typedef for a smart CRPOperator with real constant
   typedef EPMSHARED_PTR<PolynomialOperator<CRPOperator<EPMFloat> > >  SmartDRPOp;

   /// Typedef for a smart CLRPOperator with real constant
   typedef EPMSHARED_PTR<PolynomialOperator<CLRPOperator<EPMFloat> > >  SmartDLRPOp;

   /// Typedef for a smart CPOperator with complex constant
   typedef EPMSHARED_PTR<PolynomialOperator<CPOperator<EPMComplex> > >  SmartZPOp;

   /// Typedef for a smart CLPOperator with complex constant
   typedef EPMSHARED_PTR<PolynomialOperator<CLPOperator<EPMComplex> > >  SmartZLPOp;

   /// Typedef for a smart CRPOperator with complex constant
   typedef EPMSHARED_PTR<PolynomialOperator<CRPOperator<EPMComplex> > >  SmartZRPOp;

   /// Typedef for a smart CLRPOperator with complex constant
   typedef EPMSHARED_PTR<PolynomialOperator<CLRPOperator<EPMComplex> > >  SmartZLRPOp;

}

#endif // POLYOPSTYPES_HPP
