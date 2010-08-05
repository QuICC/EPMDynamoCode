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

namespace EPMDynamo {

   /**
    * @brief Some tools to create the CSCS output file
    *
    * \epmTodo The field computation need correction for the pole and origin values
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
          * @param reversed Reversed order array
          * @param data Input data
          */
         static void reverseOrder(Array &reverse, const Array &data);

         /**
          * @brief Adapt RTP field to CSCS output
          *
          * @param corrected Corrected RTP field
          * @param field Input field
          * @param totNth Total number of theta points
          * @param th0 Array of start index of theta points
          * @param nTh Array of number of theta points
          */
         static void padRTPField(std::vector<Matrix> &corrected, const std::vector<Matrix> &field, const int totNth, const ArrayI &th0, const ArrayI &nTh);

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
