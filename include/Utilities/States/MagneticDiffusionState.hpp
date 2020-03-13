/** \file MagneticDiffusionState.hpp
 *  \brief Initial to compute the magnetic diffusion rates. Initialises spherical bessel basic states
 */

#ifndef MAGNETICDIFFUSIONSTATE_HPP
#define MAGNETICDIFFUSIONSTATE_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits of Magnetic diffusion spherical Bessel state generator
    */
   class MagneticDiffusionStateTraits
   {
      public:
         /// Requires RTP Codensity computations
         static const bool UseRTPCodensity = false;

         /// Requires spectral Codensity computations
         static const bool UseSpecCodensity = false;

         /// Requires Codensity gradient computations
         static const bool UseCodensityGrad = false;

         /// Requires RTP Magnetic computations
         static const bool UseRTPMagnetic = false;

         /// Requires spectral Magnetic computations
         static const bool UseSpecMagnetic = true;

         /// Requires Magnetic curl computations
         static const bool UseMagneticCurl = false;

         /// Requires RTP Velocity computations
         static const bool UseRTPVelocity = false;

         /// Requires Velocity computations
         static const bool UseSpecVelocity = false;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;
   };

   /**
    * \brief Initial to compute the magnetic diffusion rates. Initialises spherical bessel basic states
    *
    * \epmBug State generation is not aware of parallelisation
    */
   template <typename TGenTraits> class MagneticDiffusionState
   {
      public:
         /// Typdef for the StateTraits type
         typedef MagneticDiffusionStateTraits  StateTraits;

         /**
          * @brief Initial amplitude
          */
         static const EPMFloat INITIAL_AMPLITUDE;

         /**
          * @brief Highest harmonic degree with spherical Bessel projection available
          */
         static const int SPHBESSEL_MAXL;

         /**
          * @brief Highest radial basis projection with spherical Bessel projection available
          */
         static const int SPHBESSEL_MAXN;

         /**
          * @brief Order of the precomputed basic state
          */
         static const int SPHBESSEL_ORDER;

         /**
          * @brief Initial state file extension
          */
         static const std::string  FILE_EXTENSION;

         /**
          * @brief Toroidal initial state file basename
          */
         static const std::string  TOROIDAL_FILE_BASENAME;

         /**
          * @brief Poloidal initial state file basename
          */
         static const std::string  POLOIDAL_FILE_BASENAME;

         /// Typedef for the codensity type
         typedef typename TGenTraits::CodType  Codensity;

         /// Typedef for the magnetic type
         typedef typename TGenTraits::MagType  Magnetic;

         /// Typedef for the velocity type
         typedef typename TGenTraits::VelType  Velocity;

         /**
          * @brief Set RTP value for the codensity scalar
          */
         static void setRTPCodensity(Codensity &codC);

         /**
          * @brief Set RTP value for the magnetic field
          */
         static void setRTPMagnetic(Magnetic &magB);

         /**
          * @brief Set RTP value for the velocity field
          */
         static void setRTPVelocity(Velocity &velV);

         /**
          * @brief Set spectral coefficients for the codensity scalar
          */
         static void setSpecCodensity(Codensity &codC);

         /**
          * @brief Set spectral coefficients for the velocity field
          */
         static void setSpecMagnetic(Magnetic &magB);

         /**
          * @brief Set spectral coefficients for the magnetic field
          */
         static void setSpecVelocity(Velocity &velV);

         /**
          * @brief Get the initial state for input file
          */
         static Array getPrecomputedSpectrum(const std::string &basename, const int l, const int n);

      private:
         /**
          * @brief Private constructor
          */
         MagneticDiffusionState();

         /**
          * @brief Destructor
          */
         virtual ~MagneticDiffusionState() {};
   };

   template <typename TGenTraits> const EPMFloat MagneticDiffusionState<TGenTraits>::INITIAL_AMPLITUDE = 1.0e0;

   template <typename TGenTraits> const int MagneticDiffusionState<TGenTraits>::SPHBESSEL_MAXL = 15;

   template <typename TGenTraits> const int MagneticDiffusionState<TGenTraits>::SPHBESSEL_MAXN = 15;

   template <typename TGenTraits> const int MagneticDiffusionState<TGenTraits>::SPHBESSEL_ORDER = 1;

   template <typename TGenTraits> const std::string MagneticDiffusionState<TGenTraits>::FILE_EXTENSION = ".dat";

   template <typename TGenTraits> const std::string MagneticDiffusionState<TGenTraits>::TOROIDAL_FILE_BASENAME = "ToroidalSphericalBesselSpectrum";

   template <typename TGenTraits> const std::string MagneticDiffusionState<TGenTraits>::POLOIDAL_FILE_BASENAME = "PoloidalSphericalBesselSpectrum";

   template <typename TGenTraits> void MagneticDiffusionState<TGenTraits>::setRTPCodensity(typename MagneticDiffusionState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         codC.rOc().rRTP().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void MagneticDiffusionState<TGenTraits>::setRTPMagnetic(typename MagneticDiffusionState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
         magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
      }
   }

   template <typename TGenTraits> void MagneticDiffusionState<TGenTraits>::setRTPVelocity(typename MagneticDiffusionState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
      {
         velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
         velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
      }
   }

   template <typename TGenTraits> void MagneticDiffusionState<TGenTraits>::setSpecCodensity(typename MagneticDiffusionState<TGenTraits>::Codensity &codC)
   {
      SmartTruncation pTrunc = codC.oc().trunc();

      // Set some perturbation random energy
      for(int l=0; l < pTrunc->local()->spec()->nL(); ++l)
      {
         codC.rOc().rPerturbation().rLShell(l).setConstant(EPMComplex(0.0,0.0));
      }
   }

   template <typename TGenTraits> void MagneticDiffusionState<TGenTraits>::setSpecMagnetic(typename MagneticDiffusionState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      int maxNn;
      Array init;
      int minNl = std::min(pTrunc->local()->spec()->nL(),MagneticDiffusionState<TGenTraits>::SPHBESSEL_MAXL+1);

      for(int l=1; l < minNl; ++l)
      {
      //   
      // Set some initial spherical Bessel eigen state in Toroidal component
      //
         
         // Get the precomputed state (if avaiable)
         init = MagneticDiffusionState<TGenTraits>::getPrecomputedSpectrum(MagneticDiffusionState<TGenTraits>::TOROIDAL_FILE_BASENAME, l, MagneticDiffusionState<TGenTraits>::SPHBESSEL_ORDER);

         // Check that we actually got something
         if(init.size() != 0)
         {
            maxNn = std::max(static_cast<int>(init.size()), pTrunc->sim()->rad()->nN());
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0).head(maxNn).real() = init.head(maxNn);
         }

      //
      // Set some initial spherical Bessel eigen state in Poloidal component
      //
         
         // Get the precomputed state (if avaiable)
         init = MagneticDiffusionState<TGenTraits>::getPrecomputedSpectrum(MagneticDiffusionState<TGenTraits>::POLOIDAL_FILE_BASENAME, l, MagneticDiffusionState<TGenTraits>::SPHBESSEL_ORDER);

         // Check that we actually got something
         if(init.size() != 0)
         {
            maxNn = std::max(static_cast<int>(init.size()), pTrunc->sim()->rad()->nN());
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0).head(maxNn).real() = init.head(maxNn);
         }

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag(0.0);
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag(0.0);
         }
      }
   }

   template <typename TGenTraits> void MagneticDiffusionState<TGenTraits>::setSpecVelocity(typename MagneticDiffusionState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL(); ++l)
      {
         // Set some perturbation random energy in Toroidal component
         velV.rOc().rPerturbation().rTor().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Set some perturbation random energy in Poloidal component
         velV.rOc().rPerturbation().rPol().rLShell(l).setConstant(EPMComplex(0.0,0.0));

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            velV.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag(0.0);
            velV.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag(0.0);
         }
      }
   }

   template <typename TGenTraits> MagneticDiffusionState<TGenTraits>::MagneticDiffusionState()
   {
   }

   template <typename TGenTraits> Array MagneticDiffusionState<TGenTraits>::getPrecomputedSpectrum(const std::string &basename, const int l, const int n)
   {
      std::ifstream file;

      // Create a stringstream
      std::stringstream ss;

      // Convert the requested harmonic degree into a string
      ss << l;

      // Create full filename
      std::string filename = basename + "L" + ss.str() + MagneticDiffusionState<TGenTraits>::FILE_EXTENSION;

      // Open the input file (if possible)
      file.open(filename.c_str());

      // Check for success
      Array state;
      if(file.good())
      {
         // Number of radial orders that have been stored
         int maxN;
         //Number of modes that have been stored
         int modes;

         file >> maxN;
         file >> modes;

         if(n <= modes)
         {
            // Resize state array to be able to store spectrum
            state.resize(maxN+1);

            // For simplificy simple read through the one not requested
            for(int j=1; j < n; j++)
            {
               for(int i=0; i < maxN + 1; i++)
               {
                  file >> state(0);
               }
            }

            // read in the requested state spectrum
            for(int i=0; i < maxN + 1; i++)
            {
               file >> state(i);
            }
         } else
         {
            std::cerr << "Requested modes has not been precomputed" << std::endl;
         }

         // Close the input file
         file.close();

      // Opening the file failed
      } else
      {
         std::cerr << "Couldn't open base state file: " + filename << std::endl;
      }

      // Return state array (which might be empty)
      return state;
   }

}


#endif // MAGNETICDIFFUSIONSTATE_HPP
