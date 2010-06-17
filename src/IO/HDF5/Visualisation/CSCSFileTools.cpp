/** \file CSCSFileTools.cpp
 *  \brief Source of the tools used for creating the CSCS visualisation output format
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/Visualisation/CSCSFileTools.hpp"

// Project includes
//
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   void  CSCSFileTools::padPhiGrid(Array &corrected, const Array &phi)
   {
      // Create correct size
      corrected.resize(phi.size()+1);

      // Put grid back in place and add 2*pi
      corrected.segment(0,phi.size()) = phi;
      corrected(corrected.size()-1) = 2.0*MathConstants::PI;
   }

   void CSCSFileTools::padRTPField(std::vector<Matrix> &corrected, const std::vector<Matrix> &field, const int totNth, const ArrayI &th0, const ArrayI &nTh)
   {
      // Make sure corrected is empty
      if(corrected.size() != 0)
      {
         corrected.clear();
      }

      // Simulation Phi size
      int rows = field.at(0).rows();

      // Padded Phi size
      int crows = rows + 1;

      // Padded Theta size
      int cols;

      for(unsigned int n=0; n < field.size(); ++n)
      {
         cols = field.at(n).cols();

         // Create matrix
         corrected.push_back(Matrix(crows, cols));

         // Put field values in
         corrected.at(n).block(0,0,rows,cols) = field.at(n);

         // Copy over the phi=pi values from phi=0
         corrected.at(n).row(crows-1) = corrected.at(n).row(0);

         // Reorder the theta columns
         for(int j=0; j < cols/2; ++j)
         {
            corrected.at(n).col(j).swap(corrected.at(n).col(cols - 1 - j));
         }
      }
   }

   CSCSFileTools::CSCSFileTools()
   {
   }
}
