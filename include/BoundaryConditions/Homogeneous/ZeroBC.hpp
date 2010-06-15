/** \file ZeroBC.hpp
 *  \brief Implementation of the basis Zero value boundary condition
 */

#ifndef ZEROBC_HPP
#define ZEROBC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/HomogeneousBC.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the zero boundary condition
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ZeroBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         ZeroBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~ZeroBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
   };

   template <typename TSimType> ZeroBC<TSimType>::ZeroBC(const typename ZeroBC<TSimType>::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      fillLHSBCValues(basis);
   }

   template <typename TSimType> void ZeroBC<TSimType>::fillLHSBCValues(const typename ZeroBC<TSimType>::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

}

#endif // ZEROBC_HPP
