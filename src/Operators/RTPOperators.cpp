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

   void RTPOperators::cross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const EPMFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Compute R component
            rCrossVal.rR().rShell(n) = left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array();
            // Compute Theta component
            rCrossVal.rTheta().rShell(n) = left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array();
            // Compute Phi component
            rCrossVal.rPhi().rShell(n) = left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Compute R component
            rCrossVal.rR().rShell(n) = coeff*(left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array());
            // Compute Theta component
            rCrossVal.rTheta().rShell(n) = coeff*(left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array());
            // Compute Phi component
            rCrossVal.rPhi().rShell(n) = coeff*(left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array());
         }
      }
   }

   void RTPOperators::addCross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const EPMFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Add to R component
            rCrossVal.rR().rShell(n) += left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array();
            // Add to Theta component
            rCrossVal.rTheta().rShell(n) += left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array();
            // Add to Phi component
            rCrossVal.rPhi().rShell(n) += left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Add to R component
            rCrossVal.rR().rShell(n) += coeff*(left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array());
            // Add to Theta component
            rCrossVal.rTheta().rShell(n) += coeff*(left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array());
            // Add to Phi component
            rCrossVal.rPhi().rShell(n) += coeff*(left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array());
         }
      }
   }

   void RTPOperators::subCross(RTPField &rCrossVal, const RTPField &left, const RTPField &right, const EPMFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Substract from R component
            rCrossVal.rR().rShell(n) -= left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array();
            // Substract from Theta component
            rCrossVal.rTheta().rShell(n) -= left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array();
            // Substract from Phi component
            rCrossVal.rPhi().rShell(n) -= left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Substract from R component
            rCrossVal.rR().rShell(n) -= coeff*(left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array());
            // Substract from Theta component
            rCrossVal.rTheta().rShell(n) -= coeff*(left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array());
            // Substract from Phi component
            rCrossVal.rPhi().rShell(n) -= coeff*(left.r().shell(n).arrray()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array());
         }
      }
   }

   void RTPOperators::dot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const EPMFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) = left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) = coeff*(left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n) + left.phi().shell(n).array()*right.phi().shell(n).array());
         }
      }
   }
         
   void RTPOperators::addDot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const EPMFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) += left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) += coeff*(left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array());
         }
      }
   }
         
   void RTPOperators::subDot(RTPScalar &rDot, const RTPField &left, const RTPField &right, const EPMFloat coeff)
   {
      int nR = left.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) -= left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) -= coeff*(left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array());
         }
      }
   }
         
   void RTPOperators::radVect(RTPField &rField, const RTPScalar &scalar, const EPMFloat coeff)
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
         
   void RTPOperators::addRadVect(RTPField &rField, const RTPScalar &scalar, const EPMFloat coeff)
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
         
   void RTPOperators::subRadVect(RTPField &rField, const RTPScalar &scalar, const EPMFloat coeff)
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
         
   void RTPOperators::radVectProj(RTPScalar &rScalar, const RTPScalar &scalar, const EPMFloat coeff)
   {
      int nR = scalar.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Set scalar product
            rScalar.rShell(n) = scalar.radGrid(n)*scalar.shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Set scalar product
            rScalar.rShell(n) = coeff*(scalar.radGrid(n)*scalar.shell(n));
         }
      }
   }
         
   void RTPOperators::addRadVectProj(RTPScalar &rScalar, const RTPScalar &scalar, const EPMFloat coeff)
   {
      int nR = scalar.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Add scalar product
            rScalar.rShell(n) += scalar.radGrid(n)*scalar.shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Add scalar product
            rScalar.rShell(n) += coeff*(scalar.radGrid(n)*scalar.shell(n));
         }
      }
   }
         
   void RTPOperators::subRadVectProj(RTPScalar &rScalar, const RTPScalar &scalar, const EPMFloat coeff)
   {
      int nR = scalar.nR();
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            // Substract scalar product
            rScalar.rShell(n) -= scalar.radGrid(n)*scalar.shell(n);
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Substract scalar product
            rScalar.rShell(n) -= coeff*(scalar.radGrid(n)*scalar.shell(n));
         }
      }
   }
         
   void RTPOperators::crossZVect(RTPField &rField, const RTPField &inField, const EPMFloat coeff)
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
         
   void RTPOperators::addCrossZVect(RTPField &rField, const RTPField &inField, const EPMFloat coeff)
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
         
   void RTPOperators::subCrossZVect(RTPField &rField, const RTPField &inField, const EPMFloat coeff)
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
