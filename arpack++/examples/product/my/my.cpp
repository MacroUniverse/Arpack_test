#include "arssym.h"
#include "blas1c.h"
#include "lapackc.h"
#include <cmath>

template<class ART>
class MatrixWithProduct {

 private:

  int m, n; // Number of rows and columns.

 public:

  int nrows() { return m; }

  int ncols() { return n; }

  virtual void MultMv(ART* v, ART* w) = 0;
  // Matrix-vector product: w = M*v.

  MatrixWithProduct(int nrows, int ncols = 0)
  // Constructor.
  {
    m = nrows;
    n = (ncols?ncols:nrows);
  } // Constructor.

}; // MatrixWithProduct

template<class ART>
class SymMatrixB: public MatrixWithProduct<ART> {

 private:

  ART  shift;
  int  decsize;

  void FactorDataDeallocate();

 public:

  void FactorOP();

  void MultMv(ART* v, ART* w);

  // void MultOPv(ART* v, ART* w);

  SymMatrixB(int nv);

  SymMatrixB(int nv, ART shiftv);

  virtual ~SymMatrixB();

}; // SymMatrixB.


template<class ART>
void SymMatrixB<ART>::MultMv(ART* u, ART* v)
{

  for (int i = 0; i < 9; ++i)
	v[i] = u[i+1];
  for (int i = 1; i < 10; ++i)
	v[i] += u[i-1];


} //  MultMv.


template<class ART>
inline SymMatrixB<ART>::SymMatrixB(int nval): MatrixWithProduct<ART>(nval)
// Constructor

{} // Constructor.


template<class ART>
inline SymMatrixB<ART>::
SymMatrixB(int nv, ART shiftv): MatrixWithProduct<ART>(nv)
// Constructor with shift.

{} // Constructor with shift.


template<class ART>
inline SymMatrixB<ART>::~SymMatrixB()
// Destructor

{} // Destructor.


template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARSymStdEig<ARFLOAT, ARMATRIX> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "std::cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARSymStdEig \n";
  std::cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << std::endl;
  switch (mode) {
  case 1:
    std::cout << "Regular mode" << std::endl << std::endl;
    break;
  case 3: 
    std::cout << "Shift and invert mode" << std::endl << std::endl;
  }

  std::cout << "Dimension of the system            : " << n              << std::endl;
  std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv          << std::endl;
  std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << std::endl;
  std::cout << "Number of iterations taken         : " << Prob.GetIter() << std::endl;
  std::cout << std::endl;

  if (Prob.EigenvaluesFound()) {

    // Printing eigenvalues.

    std::cout << "Eigenvalues:" << std::endl;
    for (i=0; i<nconv; i++) {
      std::cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << std::endl;
    }
    std::cout << std::endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1) / fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      std::cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    std::cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


template<class T>
void Test(T type)
{

  // Creating a symmetric matrix.

  SymMatrixB<T> A(100,0.0); // n = 100, shift = 0.0.

  // Defining what we need: the four eigenvectors of B nearest to 0.0.
  // A.MultOPv is the function that performs the product w <- OPv.

  ARSymStdEig<T, SymMatrixB<T> > dprob(A.ncols(), 4, &A, &SymMatrixB<T>::MultMv, 0.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // Test.


int main()
{

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

  // Solving a single precision problem with n = 100.

  Test((float)0.0);

} // main
