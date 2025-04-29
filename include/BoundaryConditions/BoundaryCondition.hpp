/** \file BoundaryCondition.hpp
 *  \brief Highest level implementation of a boundary condition
 */

#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief General description of a BoundaryCondition
    *
    * Definition of the boundary condition concept independently of the used basis and the type of boundary conditions
    */
   class BoundaryCondition
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nL Number of harmonic degrees
          * @param nX Number of BC values
          */
         BoundaryCondition(const int nL, const int nX);

         /**
          * @brief Destructor
          */
         virtual ~BoundaryCondition() {};

         /**
          * @brief Get boundary condition value for position i of LHS operator
          *
          * @param i Position of the value
          * @param l Harmonic degreee l related to operator op
          */
         EPMFloat getLHSBC(const int i, const int l) const;

         /**
          * @brief Get boundary condition value array of LHS operator
          *
          * @param l Harmonic degreee l related to operator op
          */
         const Array& getLHSBC(const int l) const;

         /**
          * @brief Get boundary condition value RHS, real part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         virtual EPMFloat getRealRHSBC(const int l, const int m) const = 0;

         /**
          * @brief Get boundary condition value RHS, imaginary part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         virtual EPMFloat getImagRHSBC(const int l, const int m) const = 0;

         /**
          * @brief Get boundary condition value RHS, complex value
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         virtual EPMComplex getRHSBC(const int l, const int m) const = 0;

         /**
          * @brief Is a homogeneous boundary condition
          */
         bool isHomogeneous() const;
         
      protected:
         /**
          * @brief Get number of harmonic degrees of boundary condition
          */
         int nL() const;

         /**
          * @brief Get number of coefficients in boundary condition
          */
         int nX() const;

         /**
          * @brief Set the LHSBCValues array for the given harmonic degree
          *
          * @param l Harmonic degree l
          */
         Array&   rLHSBCValues(const int l);

         /**
          * @brief Set the Predictor BC to one constant
          *
          * @param c Constant for the BC
          */
         void   setConstPredictorBC(const EPMFloat c);

         /**
          * @brief Set the Predictor BC values
          */
         Array&   rPredictorBC();

         /**
          * @brief Flag for homogeneous boundary condition
          */
         bool mIsHomogeneous;

      private:

         /**
          * @brief Number of harmonic degrees
          */
         int mNl;

         /**
          * @brief Number of BC values
          */
         int mNx;

         /**
          * @brief Storage for the LHS boundary values
          */
         std::vector<Array>   mLHSBCValues;

         /**
          * @brief Storage for the RHS boundary values
          */
         std::vector<ArrayZ>   mRHSBCValues;

         /**
          * @brief Initialise the boundaries vector of arrays
          */
         void initLHSBCValues();
   };

   inline bool BoundaryCondition::isHomogeneous() const
   {
      return this->mIsHomogeneous;
   }

   inline int BoundaryCondition::nL() const
   {
      return this->mNl;
   }

   inline int BoundaryCondition::nX() const
   {
      return this->mNx;
   }

   /// Typedef for a smart shared boundary condition
   typedef EPMSHARED_PTR<BoundaryCondition>  SmartBC;
}

#endif // BOUNDARYCONDITION_HPP
