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

   void CSCSFileTools::padRTPField(std::vector<Matrix> &corrected, const std::vector<Matrix> &field, SmartTruncation pTrunc)
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

   void CSCSFileTools::correctRadialPole(std::vector<Matrix> &corrected, const std::vector<Matrix> &r, const std::vector<Matrix> &theta, const std::vector<Matrix> &phi, SmartTruncation pTrunc)
   {
      // Get number of phi points
      int nPh = r.at(0).rows();

      int cols;

      // Loop over all radial shells
      for(unsigned int n = 0; n < corrected.size(); ++n)
      {
         // Get the number of columns
         cols = corrected.at(n).cols();

         // Correct north pole values: Set everything but phi = 0 to 0
         if(pTrunc->local()->rtp()->th0(n) == 0)
         {
 //           corrected.at(n).col(0).segment(1, nPh).setZero();
         }

         // Correct south pole values: Set everything but phi = 0 to 0
         if(pTrunc->local()->rtp()->th0(n) + pTrunc->local()->rtp()->nTh(n) == pTrunc->sim()->hoz()->nTh())
         {
//            corrected.at(n).col(cols-1).segment(1, nPh).setZero();
         }
      }

      // Correct orgin value if applicable
      if(pTrunc->local()->rtp()->r0() == 0)
      {
         cols = corrected.at(0).cols();
         if(pTrunc->local()->rtp()->th0(0) == 0)
         {
//            corrected.at(0).block(0,1, nPh+1, cols-1).setZero();
         } else
         {
//            corrected.at(0).setZero();
         }
      }
   }

   void CSCSFileTools::correctThetaPole(std::vector<Matrix> &corrected, const std::vector<Matrix> &r, const std::vector<Matrix> &theta, const std::vector<Matrix> &phi, SmartTruncation pTrunc)
   {
      // Get number of phi points
      int nPh = r.at(0).rows();

      int cols; 

      // Loop over all radial shells
      for(unsigned int n = 0; n < corrected.size(); ++n)
      {
         // Get the number of columns
         cols = corrected.at(n).cols();

         // Correct north pole values: Set everything but phi = 0 to 0
         if(pTrunc->local()->rtp()->th0(n) == 0)
         {
//            corrected.at(n).col(0).segment(1, nPh).setZero();
         }

         // Correct south pole values: Set everything but phi = 0 to 0
         if(pTrunc->local()->rtp()->th0(n) + pTrunc->local()->rtp()->nTh(n) == pTrunc->sim()->hoz()->nTh())
         {
//            corrected.at(n).col(cols-1).segment(1, nPh).setZero();
         }
      }

      // Correct orgin value if applicable
      if(pTrunc->local()->rtp()->r0() == 0)
      {
         cols = corrected.at(0).cols();
         if(pTrunc->local()->rtp()->th0(0) == 0)
         {
//            corrected.at(0).block(0,1, nPh+1, cols-1).setZero();
         } else
         {
//            corrected.at(0).setZero();
         }
      }
   }

   void CSCSFileTools::correctPhiPole(std::vector<Matrix> &corrected, const std::vector<Matrix> &r, const std::vector<Matrix> &theta, const std::vector<Matrix> &phi, SmartTruncation pTrunc)
   {
      // Get number of phi points
      int nPh = r.at(0).rows();

      int cols; 

      // Loop over all radial shells
      for(unsigned int n = 0; n < corrected.size(); ++n)
      {
         // Get the number of columns
         cols = corrected.at(n).cols();

         // Correct north pole values: Set everything but phi = 0 to 0
         if(pTrunc->local()->rtp()->th0(n) == 0)
         {
//            corrected.at(n).col(0).segment(1, nPh).setZero();
         }

         // Correct south pole values: Set everything but phi = 0 to 0
         if(pTrunc->local()->rtp()->th0(n) + pTrunc->local()->rtp()->nTh(n) == pTrunc->sim()->hoz()->nTh())
         {
//            corrected.at(n).col(cols-1).segment(1, nPh).setZero();
         }
      }

      // Correct orgin value if applicable
      if(pTrunc->local()->rtp()->r0() == 0)
      {
         cols = corrected.at(0).cols();
         if(pTrunc->local()->rtp()->th0(0) == 0)
         {
//            corrected.at(0).block(0,1, nPh+1, cols-1).setZero();
         } else
         {
//            corrected.at(0).setZero();
         }
      }
   }

   void CSCSFileTools::reverseOrder(Array& reverse, const Array &data)
   {
      // Copy data into reverse
      reverse = data;

      // Reverse order
      for(int j=0; j < reverse.rows()/2; ++j)
      {
         reverse.row(j).swap(reverse.row(reverse.rows() - 1 - j));
      }
   }

   CSCSFileTools::CSCSFileTools()
   {
   }
}
