/** \file RTPOperators.hpp
 *  \brief Implementation of some RTP space computations
 */

#ifndef RTPOPERATORS_HPP
#define RTPOPERATORS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {
  
   // Forward declaration of RTPScalar
   class RTPScalar;
   // Forward declaration of RTPField
   class RTPField;

   /**
    * \brief Implementation of some RTP space computations
    *
    * This class provides some of the often used RTP space operations: cross product, dot product, radial projection, \f$\hat{z}\f$ cross product
    *
    * \epmBug Review these computations and write thoroughful tests
    * \epmBug Replace multiple methods with a templated one using enum
    */
   class RTPOperators
   {
      public:
         /**
          * @brief Compute the cross product of the two given vector fields
          *
          * @param rCrossVal Result of cross product
          * @param left lefthand side operand
          * @param right righthand side operand
          * @param coeff Multiplicative coefficient
          */
         static void cross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const EPMFloat coeff = 1.0);

         /**
          * @brief Compute and add the cross product of the two given vector fields
          *
          * @param rCrossVal Result of cross product
          * @param left lefthand side operand
          * @param right righthand side operand
          * @param coeff Multiplicative coefficient
          */
         static void addCross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const EPMFloat coeff = 1.0);

         /**
          * @brief Compute and substract the cross product of the two given vector fields
          *
          * @param rCrossVal Result of cross product
          * @param left lefthand side operand
          * @param right righthand side operand
          * @param coeff Multiplicative coefficient
          */
         static void subCross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const EPMFloat coeff = 1.0);

         /**
          * @brief Compute dot product of the two given vector fields
          *
          * @param rDot Result of dot product
          * @param left Left-hand side operand
          * @param right Right-hand side operand
          * @param coeff Multiplicative coefficient
          */
         static void dot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const EPMFloat coeff = 1.0);

         /**
          * @brief Compute and add dot product of the two given vector fields
          *
          * @param rDot Result of dot product
          * @param left Left-hand side operand
          * @param right Right-hand side operand
          * @param coeff Multiplicative coefficient
          */
         static void addDot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const EPMFloat coeff = 1.0);

         /**
          * @brief Compute and substract dot product of the two given vector fields
          *
          * @param rDot Result of dot product
          * @param left Left-hand side operand
          * @param right Right-hand side operand
          * @param coeff Multiplicative coefficient
          */
         static void subDot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const EPMFloat coeff = 1.0);

         /**
          * @brief Multiply scalar by radial vector
          *
          * @param rField Resulting vector field
          * @param scalar Input scalar
          * @param coeff Multiplicative coefficient
          */
         static void radVect(RTPField &rField, const RTPScalar &scalar, const EPMFloat coeff = 1.0);

         /**
          * @brief Multiply scalar by radial vector and add
          *
          * @param rField Resulting vector field
          * @param scalar Input scalar
          * @param coeff Multiplicative coefficient
          */
         static void addRadVect(RTPField &rField, const RTPScalar &scalar, const EPMFloat coeff = 1.0);

         /**
          * @brief Multiply scalar by radial vector and substract
          *
          * @param rField Resulting vector field
          * @param scalar Input scalar
          * @param coeff Multiplicative coefficient
          */
         static void subRadVect(RTPField &rField, const RTPScalar &scalar, const EPMFloat coeff = 1.0);

         /**
          * @brief Scalar product with radial vector
          *
          * @param rScalar Resulting scalar field
          * @param scalar Input scalar
          * @param coeff Multiplicative coefficient
          */
         static void radVectProj(RTPScalar &rScalar, const RTPScalar &scalar, const EPMFloat coeff = 1.0);

         /**
          * @brief Scalar product with radial vector and add
          *
          * @param rScalar Resulting scalar field
          * @param scalar Input scalar
          * @param coeff Multiplicative coefficient
          */
         static void addRadVectProj(RTPScalar &rScalar, const RTPScalar &scalar, const EPMFloat coeff = 1.0);

         /**
          * @brief Scalar product with radial vector and substract
          *
          * @param rScalar Resulting scalar field
          * @param scalar Input scalar
          * @param coeff Multiplicative coefficient
          */
         static void subRadVectProj(RTPScalar &rScalar, const RTPScalar &scalar, const EPMFloat coeff = 1.0);

         /**
          * @brief Cross product of given field with \f$\hat{z}\f$ unit vector
          *
          * @param rField Resulting vector field
          * @param inField Input field
          * @param coeff Multiplicative coefficient
          */
         static void crossZVect(RTPField &rField, const RTPField &inField, const EPMFloat coeff = 1.0);

         /**
          * @brief Add cross product of given field with \f$\hat{z}\f$ unit vector
          *
          * @param rField Resulting vector field
          * @param inField Input field
          * @param coeff Multiplicative coefficient
          */
         static void addCrossZVect(RTPField &rField, const RTPField &inField, const EPMFloat coeff = 1.0);

         /**
          * @brief Substract product of given field with \f$\hat{z}\f$ unit vector
          *
          * @param rField Resulting vector field
          * @param inField Input field
          * @param coeff Multiplicative coefficient
          */
         static void subCrossZVect(RTPField &rField, const RTPField &inField, const EPMFloat coeff = 1.0);

      protected:

      private:

         /**
         * @brief Constructor
         */
         RTPOperators();

         /**
         * @brief Destructor
         */
         virtual ~RTPOperators() {};

   };

}

#endif // RTPOPERATORS_HPP
