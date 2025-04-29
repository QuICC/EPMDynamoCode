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
            rCrossVal.rR().rShell(n) = (left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array()).matrix();
            // Compute Theta component
            rCrossVal.rTheta().rShell(n) = (left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array()).matrix();
            // Compute Phi component
            rCrossVal.rPhi().rShell(n) = (left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array()).matrix();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Compute R component
            rCrossVal.rR().rShell(n) = coeff*(left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array()).matrix();
            // Compute Theta component
            rCrossVal.rTheta().rShell(n) = coeff*(left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array()).matrix();
            // Compute Phi component
            rCrossVal.rPhi().rShell(n) = coeff*(left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array()).matrix();
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
            rCrossVal.rR().rShell(n) += (left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array()).matrix();
            // Add to Theta component
            rCrossVal.rTheta().rShell(n) += (left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array()).matrix();
            // Add to Phi component
            rCrossVal.rPhi().rShell(n) += (left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array()).matrix();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Add to R component
            rCrossVal.rR().rShell(n) += coeff*(left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array()).matrix();
            // Add to Theta component
            rCrossVal.rTheta().rShell(n) += coeff*(left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array()).matrix();
            // Add to Phi component
            rCrossVal.rPhi().rShell(n) += coeff*(left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array()).matrix();
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
            rCrossVal.rR().rShell(n) -= (left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array()).matrix();
            // Substract from Theta component
            rCrossVal.rTheta().rShell(n) -= (left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array()).matrix();
            // Substract from Phi component
            rCrossVal.rPhi().rShell(n) -= (left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array()).matrix();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            // Substract from R component
            rCrossVal.rR().rShell(n) -= coeff*(left.theta().shell(n).array()*right.phi().shell(n).array() - left.phi().shell(n).array()*right.theta().shell(n).array()).matrix();
            // Substract from Theta component
            rCrossVal.rTheta().rShell(n) -= coeff*(left.phi().shell(n).array()*right.r().shell(n).array() - left.r().shell(n).array()*right.phi().shell(n).array()).matrix();
            // Substract from Phi component
            rCrossVal.rPhi().rShell(n) -= coeff*(left.r().shell(n).array()*right.theta().shell(n).array() - left.theta().shell(n).array()*right.r().shell(n).array()).matrix();
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
            rDot.rShell(n) = (left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array()).matrix();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) = coeff*(left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array()).matrix();
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
            rDot.rShell(n) += (left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array()).matrix();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) += coeff*(left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array()).matrix();
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
            rDot.rShell(n) -= (left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array()).matrix();
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            rDot.rShell(n) -= coeff*(left.r().shell(n).array()*right.r().shell(n).array() + left.theta().shell(n).array()*right.theta().shell(n).array() + left.phi().shell(n).array()*right.phi().shell(n).array()).matrix();
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
         
   void RTPOperators::precessionFrame(RTPField &rField, const RTPField &inField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat coeff)
   {
      Array sph = inField.trunc()->sim()->hoz()->phGrid().array().sin();
      Array cph = inField.trunc()->sim()->hoz()->phGrid().array().cos();

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
               rField.rR().rShell(n).col(j) = -omega*cosalpha*inField.sTh(j, n)*inField.phi().shell(n).col(j) + omega*sinalpha*(sph.array()*inField.theta().shell(n).col(j).array()).matrix() + omega*sinalpha*(inField.cTh(j, n)*cph.array()*inField.phi().shell(n).col(j).array()).matrix();
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = -omega*cosalpha*inField.cTh(j, n)*inField.phi().shell(n).col(j) - omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.phi().shell(n).col(j).array()).matrix() - omega*sinalpha*(sph.array()*inField.r().shell(n).col(j).array()).matrix();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = omega*cosalpha*inField.sTh(j, n)*inField.r().shell(n).col(j) + omega*cosalpha*inField.cTh(j, n)*inField.theta().shell(n).col(j) - omega*sinalpha*inField.cTh(j, n)*(cph.array()*inField.r().shell(n).col(j).array()).matrix() + omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.theta().shell(n).col(j).array()).matrix();
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
               rField.rR().rShell(n).col(j) = coeff*(-omega*cosalpha*inField.sTh(j, n)*inField.phi().shell(n).col(j) + omega*sinalpha*(sph.array()*inField.theta().shell(n).col(j).array()).matrix() + omega*sinalpha*(inField.cTh(j, n)*cph.array()*inField.phi().shell(n).col(j).array()).matrix());
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = coeff*(-omega*cosalpha*inField.cTh(j, n)*inField.phi().shell(n).col(j) - omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.phi().shell(n).col(j).array()).matrix() - omega*sinalpha*(sph.array()*inField.r().shell(n).col(j).array()).matrix());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = coeff*(omega*cosalpha*inField.sTh(j, n)*inField.r().shell(n).col(j) + omega*cosalpha*inField.cTh(j, n)*inField.theta().shell(n).col(j) - omega*sinalpha*inField.cTh(j, n)*(cph.array()*inField.r().shell(n).col(j).array()).matrix() + omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.theta().shell(n).col(j).array()).matrix());
            }
         }
      }
   }
         
   void RTPOperators::addPrecessionFrame(RTPField &rField, const RTPField &inField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat coeff)
   {
      Array sph = inField.trunc()->sim()->hoz()->phGrid().array().sin();
      Array cph = inField.trunc()->sim()->hoz()->phGrid().array().cos();

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
               rField.rR().rShell(n).col(j) += -omega*cosalpha*inField.sTh(j, n)*inField.phi().shell(n).col(j) + omega*sinalpha*(sph.array()*inField.theta().shell(n).col(j).array()).matrix() + omega*sinalpha*(inField.cTh(j, n)*cph.array()*inField.phi().shell(n).col(j).array()).matrix();
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) += -omega*cosalpha*inField.cTh(j, n)*inField.phi().shell(n).col(j) - omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.phi().shell(n).col(j).array()).matrix() - omega*sinalpha*(sph.array()*inField.r().shell(n).col(j).array()).matrix();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) += omega*cosalpha*inField.sTh(j, n)*inField.r().shell(n).col(j) + omega*cosalpha*inField.cTh(j, n)*inField.theta().shell(n).col(j) - omega*sinalpha*inField.cTh(j, n)*(cph.array()*inField.r().shell(n).col(j).array()).matrix() + omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.theta().shell(n).col(j).array()).matrix();
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
               rField.rR().rShell(n).col(j) += coeff*(-omega*cosalpha*inField.sTh(j, n)*inField.phi().shell(n).col(j) + omega*sinalpha*(sph.array()*inField.theta().shell(n).col(j).array()).matrix() + omega*sinalpha*(inField.cTh(j, n)*cph.array()*inField.phi().shell(n).col(j).array()).matrix());
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) += coeff*(-omega*cosalpha*inField.cTh(j, n)*inField.phi().shell(n).col(j) - omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.phi().shell(n).col(j).array()).matrix() - omega*sinalpha*(sph.array()*inField.r().shell(n).col(j).array()).matrix());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) += coeff*(omega*cosalpha*inField.sTh(j, n)*inField.r().shell(n).col(j) + omega*cosalpha*inField.cTh(j, n)*inField.theta().shell(n).col(j) - omega*sinalpha*inField.cTh(j, n)*(cph.array()*inField.r().shell(n).col(j).array()).matrix() + omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.theta().shell(n).col(j).array()).matrix());
            }
         }
      }
   }
         
   void RTPOperators::subPrecessionFrame(RTPField &rField, const RTPField &inField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat coeff)
   {
      Array sph = inField.trunc()->sim()->hoz()->phGrid().array().sin();
      Array cph = inField.trunc()->sim()->hoz()->phGrid().array().cos();

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
               rField.rR().rShell(n).col(j) -= -omega*cosalpha*inField.sTh(j, n)*inField.phi().shell(n).col(j) + omega*sinalpha*(sph.array()*inField.theta().shell(n).col(j).array()).matrix() + omega*sinalpha*(inField.cTh(j, n)*cph.array()*inField.phi().shell(n).col(j).array()).matrix();
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) -= -omega*cosalpha*inField.cTh(j, n)*inField.phi().shell(n).col(j) - omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.phi().shell(n).col(j).array()).matrix() - omega*sinalpha*(sph.array()*inField.r().shell(n).col(j).array()).matrix();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) -= omega*cosalpha*inField.sTh(j, n)*inField.r().shell(n).col(j) + omega*cosalpha*inField.cTh(j, n)*inField.theta().shell(n).col(j) - omega*sinalpha*inField.cTh(j, n)*(cph.array()*inField.r().shell(n).col(j).array()).matrix() + omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.theta().shell(n).col(j).array()).matrix();
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
               rField.rR().rShell(n).col(j) -= coeff*(-omega*cosalpha*inField.sTh(j, n)*inField.phi().shell(n).col(j) + omega*sinalpha*(sph.array()*inField.theta().shell(n).col(j).array()).matrix() + omega*sinalpha*(inField.cTh(j, n)*cph.array()*inField.phi().shell(n).col(j).array()).matrix());
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) -= coeff*(-omega*cosalpha*inField.cTh(j, n)*inField.phi().shell(n).col(j) - omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.phi().shell(n).col(j).array()).matrix() - omega*sinalpha*(sph.array()*inField.r().shell(n).col(j).array()).matrix());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) -= coeff*(omega*cosalpha*inField.sTh(j, n)*inField.r().shell(n).col(j) + omega*cosalpha*inField.cTh(j, n)*inField.theta().shell(n).col(j) - omega*sinalpha*inField.cTh(j, n)*(cph.array()*inField.r().shell(n).col(j).array()).matrix() + omega*sinalpha*inField.sTh(j, n)*(cph.array()*inField.theta().shell(n).col(j).array()).matrix());
            }
         }
      }
   }
         
   void RTPOperators::precession(RTPField &rField, const RTPField &inField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat time, const EPMFloat coeff)
   {
      Array spht = (inField.trunc()->sim()->hoz()->phGrid().array() + time).sin();
      Array cpht = (inField.trunc()->sim()->hoz()->phGrid().array() + time).cos();

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
               rField.rR().rShell(n).col(j) = omega*sinalpha*spht.array()*inField.theta().shell(n).col(j).array() + (omega*sinalpha*inField.cTh(j,n)*cpht.array()-(1 + omega*cosalpha)*inField.sTh(j,n))*inField.phi().shell(n).col(j).array();
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = (-omega*sinalpha*spht.array())*inField.r().shell(n).col(j).array() - ((1 + omega*cosalpha)*inField.cTh(j, n) + omega*sinalpha*inField.sTh(j, n)*cpht.array())*inField.phi().shell(n).col(j).array();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = ((1 + omega*cosalpha)*inField.sTh(j,n) - omega*sinalpha*inField.cTh(j,n)*cpht.array())*inField.r().shell(n).col(j).array() + ((1 + omega*cosalpha)*inField.cTh(j,n) + omega*sinalpha*inField.sTh(j,n)*cpht.array())*inField.theta().shell(n).col(j).array();
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
               rField.rR().rShell(n).col(j) = coeff*(omega*sinalpha*spht.array()*inField.theta().shell(n).col(j).array() + (omega*sinalpha*inField.cTh(j,n)*cpht.array()-(1 + omega*cosalpha)*inField.sTh(j,n))*inField.phi().shell(n).col(j).array());
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = coeff*((-omega*sinalpha*spht.array())*inField.r().shell(n).col(j).array() - ((1 + omega*cosalpha)*inField.cTh(j, n) + omega*sinalpha*inField.sTh(j, n)*cpht.array())*inField.phi().shell(n).col(j).array());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = coeff*(((1 + omega*cosalpha)*inField.sTh(j,n) - omega*sinalpha*inField.cTh(j,n)*cpht.array())*inField.r().shell(n).col(j).array() + ((1 + omega*cosalpha)*inField.cTh(j,n) + omega*sinalpha*inField.sTh(j,n)*cpht.array())*inField.theta().shell(n).col(j).array());
            }
         }
      }
   }
         
   void RTPOperators::addPrecession(RTPField &rField, const RTPField &inField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat time, const EPMFloat coeff)
   {
      Array spht = (inField.trunc()->sim()->hoz()->phGrid().array() + time).sin();
      Array cpht = (inField.trunc()->sim()->hoz()->phGrid().array() + time).cos();

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
               rField.rR().rShell(n).col(j).array() += omega*sinalpha*spht.array()*inField.theta().shell(n).col(j).array() + (omega*sinalpha*inField.cTh(j,n)*cpht.array()-(1 + omega*cosalpha)*inField.sTh(j,n))*inField.phi().shell(n).col(j).array();
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() += (-omega*sinalpha*spht.array())*inField.r().shell(n).col(j).array() - ((1 + omega*cosalpha)*inField.cTh(j, n) + omega*sinalpha*inField.sTh(j, n)*cpht.array())*inField.phi().shell(n).col(j).array();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() += ((1 + omega*cosalpha)*inField.sTh(j,n) - omega*sinalpha*inField.cTh(j,n)*cpht.array())*inField.r().shell(n).col(j).array() + ((1 + omega*cosalpha)*inField.cTh(j,n) + omega*sinalpha*inField.sTh(j,n)*cpht.array())*inField.theta().shell(n).col(j).array();
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
               rField.rR().rShell(n).col(j).array() += coeff*(omega*sinalpha*spht.array()*inField.theta().shell(n).col(j).array() + (omega*sinalpha*inField.cTh(j,n)*cpht.array()-(1 + omega*cosalpha)*inField.sTh(j,n))*inField.phi().shell(n).col(j).array());
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() += coeff*((-omega*sinalpha*spht.array())*inField.r().shell(n).col(j).array() - ((1 + omega*cosalpha)*inField.cTh(j, n) + omega*sinalpha*inField.sTh(j, n)*cpht.array())*inField.phi().shell(n).col(j).array());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() += coeff*(((1 + omega*cosalpha)*inField.sTh(j,n) - omega*sinalpha*inField.cTh(j,n)*cpht.array())*inField.r().shell(n).col(j).array() + ((1 + omega*cosalpha)*inField.cTh(j,n) + omega*sinalpha*inField.sTh(j,n)*cpht.array())*inField.theta().shell(n).col(j).array());
            }
         }
      }
   }
         
   void RTPOperators::subPrecession(RTPField &rField, const RTPField &inField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat time, const EPMFloat coeff)
   {
      Array spht = (inField.trunc()->sim()->hoz()->phGrid().array() + time).sin();
      Array cpht = (inField.trunc()->sim()->hoz()->phGrid().array() + time).cos();

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
               rField.rR().rShell(n).col(j).array() -= omega*sinalpha*spht.array()*inField.theta().shell(n).col(j).array() + (omega*sinalpha*inField.cTh(j,n)*cpht.array()-(1 + omega*cosalpha)*inField.sTh(j,n))*inField.phi().shell(n).col(j).array();
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() -= (-omega*sinalpha*spht.array())*inField.r().shell(n).col(j).array() - ((1 + omega*cosalpha)*inField.cTh(j, n) + omega*sinalpha*inField.sTh(j, n)*cpht.array())*inField.phi().shell(n).col(j).array();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() -= ((1 + omega*cosalpha)*inField.sTh(j,n) - omega*sinalpha*inField.cTh(j,n)*cpht.array())*inField.r().shell(n).col(j).array() + ((1 + omega*cosalpha)*inField.cTh(j,n) + omega*sinalpha*inField.sTh(j,n)*cpht.array())*inField.theta().shell(n).col(j).array();
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
               rField.rR().rShell(n).col(j).array() -= coeff*(omega*sinalpha*spht.array()*inField.theta().shell(n).col(j).array() + (omega*sinalpha*inField.cTh(j,n)*cpht.array()-(1 + omega*cosalpha)*inField.sTh(j,n))*inField.phi().shell(n).col(j).array());
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() -= coeff*((-omega*sinalpha*spht.array())*inField.r().shell(n).col(j).array() - ((1 + omega*cosalpha)*inField.cTh(j, n) + omega*sinalpha*inField.sTh(j, n)*cpht.array())*inField.phi().shell(n).col(j).array());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() -= coeff*(((1 + omega*cosalpha)*inField.sTh(j,n) - omega*sinalpha*inField.cTh(j,n)*cpht.array())*inField.r().shell(n).col(j).array() + ((1 + omega*cosalpha)*inField.cTh(j,n) + omega*sinalpha*inField.sTh(j,n)*cpht.array())*inField.theta().shell(n).col(j).array());
            }
         }
      }
   }
         
   void RTPOperators::poincare(RTPField &rField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat time, const EPMFloat coeff)
   {
      Array spht = (rField.trunc()->sim()->hoz()->phGrid().array() + time).sin();
      Array cpht = (rField.trunc()->sim()->hoz()->phGrid().array() + time).cos();

      int nR = rField.nR();
      int nTh;
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = rField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j).setConstant(0.0);
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = -rField.radGrid(n)*omega*sinalpha*cpht;
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = rField.radGrid(n)*omega*sinalpha*rField.cTh(j, n)*spht;
            }
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = rField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j).setConstant(0.0);
               // Compute Theta component
               rField.rTheta().rShell(n).col(j) = -coeff*(rField.radGrid(n)*omega*sinalpha*cpht);
               // Compute Phi component
               rField.rPhi().rShell(n).col(j) = coeff*(rField.radGrid(n)*omega*sinalpha*rField.cTh(j, n)*spht);
            }
         }
      }
   }
         
   void RTPOperators::addPoincare(RTPField &rField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat time, const EPMFloat coeff)
   {
      Array spht = (rField.trunc()->sim()->hoz()->phGrid().array() + time).sin();
      Array cpht = (rField.trunc()->sim()->hoz()->phGrid().array() + time).cos();

      int nR = rField.nR();
      int nTh;
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = rField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j).array() += 0.0;
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() -= rField.radGrid(n)*omega*sinalpha*cpht.array();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() += rField.radGrid(n)*omega*sinalpha*rField.cTh(j, n)*spht.array();
            }
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = rField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j).array() += 0.0;
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() -= coeff*(rField.radGrid(n)*omega*sinalpha*cpht.array());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() += coeff*(rField.radGrid(n)*omega*sinalpha*rField.cTh(j, n)*spht.array());
            }
         }
      }
   }
         
   void RTPOperators::subPoincare(RTPField &rField, const EPMFloat omega, const EPMFloat cosalpha, const EPMFloat sinalpha, const EPMFloat time, const EPMFloat coeff)
   {
      Array spht = (rField.trunc()->sim()->hoz()->phGrid().array() + time).sin();
      Array cpht = (rField.trunc()->sim()->hoz()->phGrid().array() + time).cos();

      int nR = rField.nR();
      int nTh;
      if(coeff == 1.0)
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = rField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j).array() -= 0.0;
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() += rField.radGrid(n)*omega*sinalpha*cpht.array();
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() -= rField.radGrid(n)*omega*sinalpha*rField.cTh(j, n)*spht.array();
            }
         }
      } else
      {
         for(int n=0; n < nR; ++n)
         {
            nTh = rField.nTh(n);
            for(int j=0; j < nTh; ++j)
            {
               // Compute R component
               rField.rR().rShell(n).col(j).array() -= 0.0;
               // Compute Theta component
               rField.rTheta().rShell(n).col(j).array() += coeff*(rField.radGrid(n)*omega*sinalpha*cpht.array());
               // Compute Phi component
               rField.rPhi().rShell(n).col(j).array() -= coeff*(rField.radGrid(n)*omega*sinalpha*rField.cTh(j, n)*spht.array());
            }
         }
      }
   }
         
}
