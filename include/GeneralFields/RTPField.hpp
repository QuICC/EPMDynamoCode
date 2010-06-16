/** \file RTPField.hpp
 *  \brief Implementation of a RTP space vector field
 */

#ifndef RTPFIELD_HPP
#define RTPFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "GeneralScalars/RTPBase.hpp"
#include "GeneralScalars/RTPScalar.hpp"
#include "Operators/RTPOperators.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a RTP space vector field
    *
    * This is no more than a set of three RTPScalar which build the vector field.
    */
   class RTPField: public RTPBase
   {
      public:
         /**
         * @brief Construct the field with the given SpectralTruncation pointer
         *
         * @param pTrunc Truncation information
         */
         RTPField(const SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~RTPField() {};

         /**
          * @brief Get the radial component
          */
         const RTPScalar&  r() const;

         /**
          * @brief Set the radial component
          */
         RTPScalar&  rR();

         /**
          * @brief Get the theta component
          */
         const RTPScalar&  theta() const;

         /**
          * @brief Set the theta component
          */
         RTPScalar&  rTheta();

         /**
          * @brief Get the phi component
          */
         const RTPScalar&  phi() const;

         /**
          * @brief Set the phi component
          */
         RTPScalar&  rPhi();

         /**
          * @brief Compute the cross product
          *
          * @param rCross Result of cross product
          * @param right Right operand of dot product
          * @param coeff Multiplicative coefficient (with default value at 1.0)
          *
          * \tparam TOp Type of operation: (0) Set result, (1) Add result, (-1) Substract result
          */
         template <int TOp> void cross(RTPField &rCross, const RTPField &right, const DynamoFloat coeff = 1.0) const;

         /**
          * @brief Compute the dot product
          *
          * @param rDot Result of dot product
          * @param right Right operand of dot product
          * @param coeff Multiplicative coefficient (with default value at 1.0)
          *
          * \tparam TOp Type of operation: (0) Set result, (1) Add result, (-1) Substract result
          */
         template <int TOp> void dot(RTPScalar &rDot, const RTPField &right, const DynamoFloat coeff = 1.0) const;

         /**
          * @brief Compute the cross product with \f$\hat{z}\f$
          *
          * @param rField Resulting field
          * @param coeff Multiplicative coefficient (with default value at 1.0)
          *
          * \tparam TOp Type of operation: (0) Set result, (1) Add result, (-1) Substract result
          */
         template <int TOp> void crossZVect(RTPField &rField, const DynamoFloat coeff = 1.0) const;

      protected:

      private:
         /**
          * @brief Storage for radial component
          */
         RTPScalar   mRComponent;

         /**
          * @brief Storage for \f$\theta\f$ component
          */
         RTPScalar   mThComponent;

         /**
          * @brief Storage for \f$\phi\f$ component
          */
         RTPScalar   mPhComponent;
   };

   inline const RTPScalar& RTPField::r() const
   {
      return this->mRComponent;
   }

   inline RTPScalar& RTPField::rR()
   {
      return this->mRComponent;
   }

   inline const RTPScalar& RTPField::theta() const
   {
      return this->mThComponent;
   }

   inline RTPScalar& RTPField::rTheta()
   {
      return this->mThComponent;
   }

   inline const RTPScalar& RTPField::phi() const
   {
      return this->mPhComponent;
   }

   inline RTPScalar& RTPField::rPhi()
   {
      return this->mPhComponent;
   }

   template <int TOp> inline void RTPField::cross(RTPField &rCross, const RTPField &right, const DynamoFloat coeff) const
   {
      if(TOp == 0)
      {
         RTPOperators::cross(rCross, *this, right, coeff);
      }
      else if(TOp > 0)
      {
         RTPOperators::addCross(rCross, *this, right, coeff);
      }
      else if(TOp < 1)
      {
         RTPOperators::subCross(rCross, *this, right, coeff);
      }
   }

   template <int TOp> inline void RTPField::dot(RTPScalar &rDot, const RTPField &right, const DynamoFloat coeff) const
   {
      if(TOp == 0)
      {
         RTPOperators::dot(rDot, *this, right, coeff);
      }
      else if(TOp > 0)
      {
         RTPOperators::addDot(rDot, *this, right, coeff);
      }
      else if(TOp < 1)
      {
         RTPOperators::subDot(rDot, *this, right, coeff);
      }
   }

   template <int TOp> inline void RTPField::crossZVect(RTPField &rField, const DynamoFloat coeff) const
   {
      if(TOp == 0)
      {
         RTPOperators::crossZVect(rField, *this, coeff);
      }
      else if(TOp > 0)
      {
         RTPOperators::addCrossZVect(rField, *this, coeff);
      }
      else if(TOp < 0)
      {
         RTPOperators::subCrossZVect(rField, *this, coeff);
      }
   }
}

#endif // RTPFIELD_HPP
