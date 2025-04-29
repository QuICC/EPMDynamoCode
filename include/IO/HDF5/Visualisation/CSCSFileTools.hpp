/** \file CSCSFileTools.hpp
 *  \brief Implementation of some tools required to create de CSCS visualisation output format
 */

#ifndef CSCSFILETOOLS_HPP
#define CSCSFILETOOLS_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Some tools to create the CSCS output file
    *
    * \epmTodo the corrections for the pole and origin values could be improved
    * \epmTodo Doesn not work yet with tubular parallelisation
    */
   class CSCSFileTools
   {
      public:
         /**
          * @brief Adapt phi grid to CSCS output
          *
          * @param corrected Corrected phi grid values
          * @param phi Input phi grid
          */
         static void padPhiGrid(Array &corrected, const Array &phi);

         /**
          * @brief Reverse the order of the given array
          *
          * @param reverse Reverse order array
          * @param data Input data
          */
         static void reverseOrder(Array &reverse, const Array &data);

         /**
          * @brief Adapt RTP field to CSCS output
          *
          * @param corrected Corrected RTP field
          * @param field Input field
          * @param pTrunc Truncation information
          */
         static void padRTPField(std::vector<Matrix> &corrected, const std::vector<Matrix> &field, SmartTruncation pTrunc);

         /**
          * @brief Correct the pole values of the radial component
          *
          * @param corrected Corrected RTP field
          * @param r Radial component
          * @param theta Theta component
          * @param phi Phi component
          * @param pTrunc Truncation information
          */
         static void correctRadialPole(std::vector<Matrix> &corrected, const std::vector<Matrix> &r, const std::vector<Matrix> &theta, const std::vector<Matrix> &phi, SmartTruncation pTrunc);

         /**
          * @brief Correct the pole values of the theta component
          *
          * @param corrected Corrected RTP field
          * @param r Radial component
          * @param theta Theta component
          * @param phi Phi component
          * @param pTrunc Truncation information
          */
         static void correctThetaPole(std::vector<Matrix> &corrected, const std::vector<Matrix> &r, const std::vector<Matrix> &theta, const std::vector<Matrix> &phi, SmartTruncation pTrunc);

         /**
          * @brief Correct the pole values of the phi component
          *
          * @param corrected Corrected RTP field
          * @param r Radial component
          * @param theta Theta component
          * @param phi Phi component
          * @param pTrunc Truncation information
          */
         static void correctPhiPole(std::vector<Matrix> &corrected, const std::vector<Matrix> &r, const std::vector<Matrix> &theta, const std::vector<Matrix> &phi, SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~CSCSFileTools() {};

      private:
         /**
         * @brief Empty destructor
         */
         CSCSFileTools();
   };

}

#endif // CSCSFILETOOLS_HPP
