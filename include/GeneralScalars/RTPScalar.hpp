/** \file RTPScalar.hpp
 *  \brief Implementation of a RTP space scalar field data type
 */

#ifndef RTPSCALAR_HPP
#define RTPSCALAR_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/RTPBase.hpp"
#include "Operators/RTPOperators.hpp"

namespace EPMDynamo {

   // Forward declaration of RTPField
   class RTPField;

   /**
    * @brief Real space data abstraction for scalar field
    *
    * The scalar will be stored as a bunch of SphericalShell.
    */
   class RTPScalar : public RTPBase
   {
      public:
         /**
         * @brief Constructs the scalar based on the given truncation pointer
         *
         * @param pTrunc truncation information
         */
         RTPScalar(const SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~RTPScalar() {};

         /**
          * @brief Set the vector of SphericalShell
          */
         std::vector<SphericalShell>& data();

         /**
          * @brief Get the vector of SphericalShell
          */
         const std::vector<SphericalShell>& data() const;

         /**
          * @brief Get the spherical shell for a given radius index
          *
          * @param iR   Index of the radius coordinate
          */
         const SphericalShell&    shell(const int iR) const;

         /**
          * @brief Set the spherical shell at the given radius index
          *
          * @param iR   Index of the radius coordinate
          */
         SphericalShell&    rShell(const int iR);

         /**
          * @brief Multiply scalar field by radial vector
          *
          * @param rField  Resulting vector field
          * @param coeff Multiplicative coefficient
          */
         template <int TOp> void radVect(RTPField &rField, const EPMFloat = 1) const;

         /**
          * @brief Scalar product of scalar field with radial vector
          *
          * @param rScalar  Resulting scalar field
          * @param coeff Multiplicative coefficient
          */
         template <int TOp> void radVectProj(RTPScalar &rScalar, const EPMFloat = 1) const;

         /**
          * @brief Initialise values to zero
          */
         void initialiseZeros();
         
      protected:

      private:
         /**
          * @brief Resolution information
          */
         SmartTruncation   mpTrunc;

         /**
          * @brief Storage for the values of the scalar field
          *
          * The values are stored by shells in a std::vector. This has been chosen 
          * because most of the accesses are done shell by shell.
          */
         std::vector<SphericalShell>   mValues;

         /**
          * @brief Initialise the shells vector
          */
         void initShells();
   };

   inline std::vector<SphericalShell>& RTPScalar::data()
   {
      return this->mValues;
   }

   inline const std::vector<SphericalShell>& RTPScalar::data() const
   {
      return this->mValues;
   }

   inline const SphericalShell&   RTPScalar::shell(const int iR) const
   {
      return this->mValues.at(iR);
   }

   inline SphericalShell&   RTPScalar::rShell(const int iR)
   {
      return this->mValues.at(iR);
   }

   template <int TOp> inline void RTPScalar::radVect(RTPField &rField, const EPMFloat coeff) const
   {
      if(TOp == 0)
      {
         RTPOperators::radVect(rField, *this, coeff);
      }
      else if(TOp > 0)
      {
         RTPOperators::addRadVect(rField, *this, coeff);
      }
      else if(TOp < 0)
      {
         RTPOperators::subRadVect(rField, *this, coeff);
      }

   }

   template <int TOp> inline void RTPScalar::radVectProj(RTPScalar &rScalar, const EPMFloat coeff) const
   {
      if(TOp == 0)
      {
         RTPOperators::radVectProj(rScalar, *this, coeff);
      }
      else if(TOp > 0)
      {
         RTPOperators::addRadVectProj(rScalar, *this, coeff);
      }
      else if(TOp < 0)
      {
         RTPOperators::subRadVectProj(rScalar, *this, coeff);
      }

   }

}

#endif // RTPSCALAR_HPP
