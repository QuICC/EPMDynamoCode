/** \file LAPACK_Iface.hpp
 *  This file contains the declaration needed to use LAPACK routines in the code (until eigen has all needed things...)
 */

#ifndef LAPACK_IFACE_HPP
#define LAPACK_IFACE_HPP

extern "C"
{
   /**
    * @brief Declaration for the LAPACK dgetrf routine
    */
   void dgetrf_(int* m, int* n, double* A, int* lda, int* ipiv, int* info);

   /**
    * @brief Declaration for the LAPACK dgetri routine
    */
   void dgetri_(int* n, double* A, int* lda, int* ipiv, double* work, int* lwork, int* info);
   /**
    * @brief Declaration for the LAPACK dstevx routine
    */
   void dstevx_(char* jobz, char* range, int* N, double* d, double* e, double* vl, double* vu, int* il, int *iu, double* abstol, int* m, double* w, double* z, int* ldz, double* work, int* iwork, int* ifail, int* info); 

   /**
    * @brief Declaration for the LAPACK dlamch routine
    */
   double dlamch_(char* cmach);

   /**
    * @brief Declaration for the LAPACK dgetrs routine
    */
   void dgetrs_(char* trans, int* N, int* NRHS, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info); 

   /**
    * @brief Declaration for the LAPACK dtrtrs routine
    */
   void dtrtrs_(char* uplo, char* trans, char* diag, int* N, int* NRHS, double* A, int* lda, double* B, int* ldb, int* info); 

   /**
    * @brief Declaration for the LAPACK dgeev routine
    */
   void dgeev_(char* jobvl, char* jobvr, int* N, double* A, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info); 
}

#endif // LAPACK_IFACE_HPP
