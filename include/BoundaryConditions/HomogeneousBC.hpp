/** \file HomogeneousBC.hpp
 *  \brief General implemenation for a homogeneous boundary condition
 */

#ifndef HOMOGENEOUSBC_HPP
#define HOMOGENEOUSBC_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"

namespace EPMDynamo {

   /**
    * @brief General description of a HomogeneousBC
    *
    * This class defines the HomogeneousBC concept independently of the used basis
    */
   class HomogeneousBC: public BoundaryCondition
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nL Number of harmonic degrees
          * @param nX Number of BC values
          */
         HomogeneousBC(const int nL, const int nX);

         /**
          * @brief Destructor
          */
         virtual ~HomogeneousBC() {};

         /**
          * @brief Get boundary condition value RHS, real part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         DynamoFloat getRealRHSBC(const int l, const int m) const;

         /**
          * @brief Get boundary condition value RHS, imaginary part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         DynamoFloat getImagRHSBC(const int l, const int m) const;

         /**
          * @brief Get boundary condition value RHS, complex value
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         DynamoComplex getRHSBC(const int l, const int m) const;
         
      protected:

      private:
   };

   inline DynamoFloat HomogeneousBC::getRealRHSBC(const int l, const int m) const
   {
      return 0.0;
   }

   inline DynamoFloat HomogeneousBC::getImagRHSBC(const int l, const int m) const
   {
      return 0.0;
   }

   inline DynamoComplex HomogeneousBC::getRHSBC(const int l, const int m) const
   {
      return DynamoComplex(0.0);
   }
}

#endif // HOMOGENEOUSBC_HPP
