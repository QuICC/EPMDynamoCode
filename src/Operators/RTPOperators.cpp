/** \file RTPOperators.cpp
 *  \brief Source of the implementation of some often used RTP space operations
 */

// System includes
//

// External includes
//

// Class include
//
#include "Operators/RTPOperators.hpp"

// Project includes
//
#include "GeneralFields/RTPField.hpp"

namespace EPMDynamo {

   RTPOperators::RTPOperators()
   {
   }

   void RTPOperators::cross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const DynamoFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Compute R component
            rCrossVal.rR().rShell(n) = left.theta().shell(n).cwise()*right.phi().shell(n) - left.phi().shell(n).cwise()*right.theta().shell(n);
            // Compute Theta component
            rCrossVal.rTheta().rShell(n) = left.phi().shell(n).cwise()*right.r().shell(n) - left.r().shell(n).cwise()*right.phi().shell(n);
            // Compute Phi component
            rCrossVal.rPhi().rShell(n) = left.r().shell(n).cwise()*right.theta().shell(n) - left.theta().shell(n).cwise()*right.r().shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Compute R component
            rCrossVal.rR().rShell(n) = coeff*(left.theta().shell(n).cwise()*right.phi().shell(n) - left.phi().shell(n).cwise()*right.theta().shell(n));
            // Compute Theta component
            rCrossVal.rTheta().rShell(n) = coeff*(left.phi().shell(n).cwise()*right.r().shell(n) - left.r().shell(n).cwise()*right.phi().shell(n));
            // Compute Phi component
            rCrossVal.rPhi().rShell(n) = coeff*(left.r().shell(n).cwise()*right.theta().shell(n) - left.theta().shell(n).cwise()*right.r().shell(n));
         }
      }
   }

   void RTPOperators::addCross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const DynamoFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Add to R component
            rCrossVal.rR().rShell(n) += left.theta().shell(n).cwise()*right.phi().shell(n) - left.phi().shell(n).cwise()*right.theta().shell(n);
            // Add to Theta component
            rCrossVal.rTheta().rShell(n) += left.phi().shell(n).cwise()*right.r().shell(n) - left.r().shell(n).cwise()*right.phi().shell(n);
            // Add to Phi component
            rCrossVal.rPhi().rShell(n) += left.r().shell(n).cwise()*right.theta().shell(n) - left.theta().shell(n).cwise()*right.r().shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Add to R component
            rCrossVal.rR().rShell(n) += coeff*(left.theta().shell(n).cwise()*right.phi().shell(n) - left.phi().shell(n).cwise()*right.theta().shell(n));
            // Add to Theta component
            rCrossVal.rTheta().rShell(n) += coeff*(left.phi().shell(n).cwise()*right.r().shell(n) - left.r().shell(n).cwise()*right.phi().shell(n));
            // Add to Phi component
            rCrossVal.rPhi().rShell(n) += coeff*(left.r().shell(n).cwise()*right.theta().shell(n) - left.theta().shell(n).cwise()*right.r().shell(n));
         }
      }
   }

   void RTPOperators::subCross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const DynamoFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Substract from R component
            rCrossVal.rR().rShell(n) -= left.theta().shell(n).cwise()*right.phi().shell(n) - left.phi().shell(n).cwise()*right.theta().shell(n);
            // Substract from Theta component
            rCrossVal.rTheta().rShell(n) -= left.phi().shell(n).cwise()*right.r().shell(n) - left.r().shell(n).cwise()*right.phi().shell(n);
            // Substract from Phi component
            rCrossVal.rPhi().rShell(n) -= left.r().shell(n).cwise()*right.theta().shell(n) - left.theta().shell(n).cwise()*right.r().shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Substract from R component
            rCrossVal.rR().rShell(n) -= coeff*(left.theta().shell(n).cwise()*right.phi().shell(n) - left.phi().shell(n).cwise()*right.theta().shell(n));
            // Substract from Theta component
            rCrossVal.rTheta().rShell(n) -= coeff*(left.phi().shell(n).cwise()*right.r().shell(n) - left.r().shell(n).cwise()*right.phi().shell(n));
            // Substract from Phi component
            rCrossVal.rPhi().rShell(n) -= coeff*(left.r().shell(n).cwise()*right.theta().shell(n) - left.theta().shell(n).cwise()*right.r().shell(n));
         }
      }
   }

   void RTPOperators::dot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const DynamoFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) = left.r().shell(n).cwise()*right.r().shell(n) + left.theta().shell(n).cwise()*right.theta().shell(n) + left.phi().shell(n).cwise()*right.phi().shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) = coeff*(left.r().shell(n).cwise()*right.r().shell(n) + left.theta().shell(n).cwise()*right.theta().shell(n) + left.phi().shell(n).cwise()*right.phi().shell(n));
         }
      }
   }
         
   void RTPOperators::addDot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const DynamoFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) += left.r().shell(n).cwise()*right.r().shell(n) + left.theta().shell(n).cwise()*right.theta().shell(n) + left.phi().shell(n).cwise()*right.phi().shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) += coeff*(left.r().shell(n).cwise()*right.r().shell(n) + left.theta().shell(n).cwise()*right.theta().shell(n) + left.phi().shell(n).cwise()*right.phi().shell(n));
         }
      }
   }
         
   void RTPOperators::subDot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const DynamoFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) -= left.r().shell(n).cwise()*right.r().shell(n) + left.theta().shell(n).cwise()*right.theta().shell(n) + left.phi().shell(n).cwise()*right.phi().shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) -= coeff*(left.r().shell(n).cwise()*right.r().shell(n) + left.theta().shell(n).cwise()*right.theta().shell(n) + left.phi().shell(n).cwise()*right.phi().shell(n));
         }
      }
   }
         
   void RTPOperators::radVect(RTPField &rField, const RTPScalar &scalar, const DynamoFloat coeff)
   {
      int nR = scalar.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Set radial component
            rField.rR().rShell(n) = scalar.radGrid(n)*scalar.shell(n);
            // Set other components to zero
            rField.rTheta().rShell(n).setConstant(0.0);
            rField.rPhi().rShell(n).setConstant(0.0);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Set radial component
            rField.rR().rShell(n) = coeff*(scalar.radGrid(n)*scalar.shell(n));
            // Set other components to zero
            rField.rTheta().rShell(n).setConstant(0.0);
            rField.rPhi().rShell(n).setConstant(0.0);
         }
      }
   }
         
   void RTPOperators::addRadVect(RTPField &rField, const RTPScalar &scalar, const DynamoFloat coeff)
   {
      int nR = scalar.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Set radial component
            rField.rR().rShell(n) += scalar.radGrid(n)*scalar.shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Set radial component
            rField.rR().rShell(n) += coeff*(scalar.radGrid(n)*scalar.shell(n));
         }
      }
   }
         
   void RTPOperators::subRadVect(RTPField &rField, const RTPScalar &scalar, const DynamoFloat coeff)
   {
      int nR = scalar.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Set radial component
            rField.rR().rShell(n) -= scalar.radGrid(n)*scalar.shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Set radial component
            rField.rR().rShell(n) -= coeff*(scalar.radGrid(n)*scalar.shell(n));
         }
      }
   }
         
   void RTPOperators::crossZVect(RTPField &rField, const RTPField &inField, const DynamoFloat coeff)
   {
      int nR = inField.nR();
      int nTh;
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = inField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j) = (-inField.sTh(j,n))*inField.phi().shell(n).col(j);
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = (-inField.cTh(j, n))*inField.phi().shell(n).col(j);
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = inField.cTh(j,n)*inField.theta().shell(n).col(j) + inField.sTh(j,n)*inField.r().shell(n).col(j);
            }
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = inField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j) = coeff*((-inField.sTh(j,n))*inField.phi().shell(n).col(j));
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = coeff*((-inField.cTh(j,n))*inField.phi().shell(n).col(j));
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = coeff*((inField.cTh(j,n))*inField.theta().shell(n).col(j) + inField.sTh(j,n)*inField.r().shell(n).col(j));
            }
         }
      }
   }
         
   void RTPOperators::addCrossZVect(RTPField &rField, const RTPField &inField, const DynamoFloat coeff)
   {
      int nR = inField.nR();
      int nTh;
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = inField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j) -= inField.sTh(j,n)*inField.phi().shell(n).col(j);
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) -= inField.cTh(j,n)*inField.phi().shell(n).col(j);
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) += inField.cTh(j,n)*inField.theta().shell(n).col(j) + inField.sTh(j,n)*inField.r().shell(n).col(j);
            }
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = inField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j) -= coeff*(inField.sTh(j,n)*inField.phi().shell(n).col(j));
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) -= coeff*(inField.cTh(j,n)*inField.phi().shell(n).col(j));
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) += coeff*(inField.cTh(j,n)*inField.theta().shell(n).col(j) + inField.sTh(j,n)*inField.r().shell(n).col(j));
            }
         }
      }
   }
         
   void RTPOperators::subCrossZVect(RTPField &rField, const RTPField &inField, const DynamoFloat coeff)
   {
      int nR = inField.nR();
      int nTh;
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = inField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j) += inField.sTh(j,n)*inField.phi().shell(n).col(j);
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) += inField.cTh(j,n)*inField.phi().shell(n).col(j);
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) -= inField.cTh(j,n)*inField.theta().shell(n).col(j) + inField.sTh(j,n)*inField.r().shell(n).col(j);
            }
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = inField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j) += coeff*(inField.sTh(j,n)*inField.phi().shell(n).col(j));
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) += coeff*(inField.cTh(j,n)*inField.phi().shell(n).col(j));
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) -= coeff*(inField.cTh(j,n)*inField.theta().shell(n).col(j) + inField.sTh(j,n)*inField.r().shell(n).col(j));
            }
         }
      }
   }
}
