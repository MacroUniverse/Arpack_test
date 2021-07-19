/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLNSMat.h.
   Arpack++ class ARluNonSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arlnspen.h"

#ifndef ARLNSMAT_H
#define ARLNSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "blas1c.h"
#include "superluc.h"
#include "arlspdef.h"
#include "arlutil.h"

template<class AR_T, class AR_S> class ARluNonSymPencil;

template<class ARTYPE, class ARFLOAT>
class ARluNonSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARluNonSymPencil<ARTYPE, ARFLOAT>;
  friend class ARluNonSymPencil<ARFLOAT, ARFLOAT>;

 protected:

  bool        factored;
  int         order;
  int         nnz;
  int*        irow;
  int*        pcol;
  int*        permc;
  int*        permr;
  double      threshold;
  ARTYPE*     a;
  SuperMatrix A;
  SuperMatrix L;
  SuperMatrix U;
  ARhbMatrix<int, ARTYPE> mat;

  bool DataOK();

  virtual void Copy(const ARluNonSymMatrix& other);

  void ClearMem();

  void SubtractAsI(ARTYPE sigma, NCformat& A, NCformat& AsI);

 public:

  int nzeros() { return nnz; }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultMtv(ARTYPE* v, ARTYPE* w);

  void MultMtMv(ARTYPE* v, ARTYPE* w);

  void MultMMtv(ARTYPE* v, ARTYPE* w);

  void Mult0MMt0v(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                    int* pcolp, double thresholdp = 0.1,
                    int orderp = 1, bool check = true);   // Square matrix.

  void DefineMatrix(int mp, int np, int nnzp, ARTYPE* ap,
                    int* irowp, int* pcolp);              // Rectangular matrix.

  ARluNonSymMatrix();
  // Short constructor that does nothing.

  ARluNonSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
                   double thresholdp = 0.1, int orderp = 1, bool check = true);
  // Long constructor (square matrix).

  ARluNonSymMatrix(int mp, int np, int nnzp, ARTYPE* ap, int* irowp,int* pcolp);
  // Long constructor (rectangular matrix).

  ARluNonSymMatrix(char* name, double thresholdp = 0.1, 
                   int orderp = 1, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARluNonSymMatrix(const ARluNonSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluNonSymMatrix() { ClearMem(); }
  // Destructor.

  ARluNonSymMatrix& operator=(const ARluNonSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARluNonSymMatrix member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARTYPE, class ARFLOAT>
bool ARluNonSymMatrix<ARTYPE, ARFLOAT>::DataOK()
{

  int i, j, k;

  // Checking if pcol is in ascending order.

  i = 0;
  while ((i!=n)&&(pcol[i]<=pcol[i+1])) i++;
  if (i!=n) return false;

  // Checking if irow components are in order and within bounds.

  for (i=0; i!=n; i++) {
    j = pcol[i];
    k = pcol[i+1]-1;
    if (j<=k) {
      if ((irow[j]<0)||(irow[k]>=n)) return false;
      while ((j!=k)&&(irow[j]<irow[j+1])) j++;
      if (j!=k) return false;
    }
  }

  return true;

} // DataOK.


template<class ARTYPE, class ARFLOAT>
inline void ARluNonSymMatrix<ARTYPE, ARFLOAT>::
Copy(const ARluNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  // Copying very fundamental variables.

  defined   = other.defined;
  factored  = other.factored;

  // Returning from here if "other" was not initialized.

  if (!defined) return;

  // Copying user-defined parameters.

  if (other.n == other.m) {
    DefineMatrix(other.n, other.nnz, other.a, other.irow,
                 other.pcol, other.threshold, other.order);
  }
  else {
    DefineMatrix(other.m, other.n, other.nnz, 
                 other.a, other.irow, other.pcol);
  }

  // Throwing the original factorization away (this procedure 
  // is really awkward, but it is necessary because there
  // is no copy function for matrices L and U in the SuperLU 
  // library and it is not a good idea to do this kind of deep 
  // copy here).

  if (factored) {
    ArpackError(ArpackError::LAPACK_ERROR, "ARluNonSymMatrix");
    factored = false;
  }

} // Copy.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::ClearMem()
{

  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree();
  }
  if (defined) {
    Destroy_SuperMatrix_Store(&A); // delete A.Store;
    delete[] permc;
    delete[] permr;
    permc = NULL;
    permr = NULL;
    A.Store = NULL;
  }

} // ClearMem.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::
SubtractAsI(ARTYPE sigma, NCformat& A, NCformat& AsI)
{

  // Defining local variables.

  int     i, j, k, end;
  ARTYPE* anzval;
  ARTYPE* inzval;

  // Telling the compiler that nzval must be viewed as a vector of ARTYPE.

  anzval = (ARTYPE*)A.nzval;
  inzval = (ARTYPE*)AsI.nzval;

  // Subtracting sigma from diagonal elements.

  k = 0;
  AsI.colptr[0] = 0;

  for (i=0; i!=n; i++) {

    j = A.colptr[i];
    end = A.colptr[i+1];

    // Copying superdiagonal elements of column i.

    while ((A.rowind[j] < i)&&(j < end)) {
      inzval[k] = anzval[j];
      AsI.rowind[k++] = A.rowind[j++];
    }

    // Verifying if A(i,i) exists.

    if ((A.rowind[j] == i)&&(j < end)) { // A(i,i) exists, subtracting sigma.
      inzval[k] = anzval[j++] - sigma;
    }
    else {                               // A(i,i) does not exist.
      inzval[k] = -sigma;
    }
    AsI.rowind[k++] = i;

    // Copying subdiagonal elements of column i.

    while (j < end ) {
      inzval[k] = anzval[j];
      AsI.rowind[k++] = A.rowind[j++];
    }

    AsI.colptr[i+1] = k;

  }

  AsI.nnz = AsI.colptr[n];

} // SubtractAsI.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::FactorA()
{

  // Defining local variables.

  int         info;
  int*        etree;
  SuperMatrix AC;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluNonSymMatrix::FactorA");
  }

  // Quitting the function if A is not square.

  if (m != n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARluNonSymMatrix::FactorA");
  }

  // Deleting previous versions of L and U.
  
  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree();
  }

  // Setting default values for gstrf parameters.

  double drop_tol        = 0.0;
  int    panel_size      = sp_ienv(1);
  int    relax           = sp_ienv(2);

  // Reserving memory for etree (used in matrix decomposition).

  etree = new int[n];

  // Defining LUStat.

  StatInit(panel_size, relax);

  // Defining the column permutation of matrix A
  // (using minimum degree ordering on A'*A).

  get_perm_c(order, &A, permc);

  // Permuting columns of A and
  // creating the elimination tree of A'*A.

  sp_preorder("N", &A, permc, etree, &AC);

  // Decomposing A.

  gstrf("N",&AC, threshold, drop_tol, relax, panel_size, etree,
        NULL, 0, permr, permc, &L, &U, &info);

  // Deleting AC and etree.

  Destroy_CompCol_Permuted(&AC);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {        // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluNonSymMatrix::FactorA");
  }
  else if (info > n) {    // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluNonSymMatrix::FactorA");
  }
  else if (info > 0) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluNonSymMatrix::FactorA");
  }

} // FactorA.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARluNonSymMatrix::FactorAsI");
  }

  // Quitting the function if A is not square.

  if (m != n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARluNonSymMatrix::FactorAsI");
  }

  // Defining local variables.

  int         info;
  int*        etree;
  int*        irowi;
  int*        pcoli;
  ARTYPE*     asi;
  SuperMatrix AsI;
  SuperMatrix AC;
  NCformat*   Astore;
  NCformat*   AsIstore;

  // Deleting previous versions of L and U.
  
  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree();
  }

  // Setting default values for gstrf parameters.

  double drop_tol        = 0.0;
  int    panel_size      = sp_ienv(1);
  int    relax           = sp_ienv(2);

  // Creating a temporary matrix AsI.

  irowi = new int[nnz+n];
  pcoli = new int[n+1];
  asi   = new ARTYPE[nnz+n];
  Create_CompCol_Matrix(&AsI, n,  n, nnz, asi, irowi, pcoli, NC, GE);

  // Subtracting sigma*I from A and storing the result on AsI.

  Astore   = (NCformat*)A.Store;
  AsIstore = (NCformat*)AsI.Store;
  SubtractAsI(sigma, *Astore, *AsIstore);

  // Reserving memory for etree (used in matrix decomposition).

  etree = new int[n];

  // Defining LUStat.

  StatInit(panel_size, relax);

  // Defining the column permutation of matrix AsI
  // (using minimum degree ordering on AsI'*AsI).

  get_perm_c(order, &AsI, permc);

  // Permuting columns of AsI and
  // creating the elimination tree of AsI'*AsI.

  sp_preorder("N", &AsI, permc, etree, &AC);

  // Decomposing AsI.

  gstrf("N",&AC, threshold, drop_tol, relax, panel_size, etree,
        NULL, 0, permr, permc, &L, &U, &info);

  // Deleting AC, AsI and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&AsI);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {        // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluNonSymMatrix::FactorAsI");
  }
  else if (info > n) {    // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluNonSymMatrix::FactorAsI");
  }
  else if (info > 0) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluNonSymMatrix::FactorAsI");
  }

} // FactorAsI.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::MultMv(ARTYPE* v, ARTYPE* w)
{

  int    i,j;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluNonSymMatrix::MultMv");
  }

  // Determining w = M.v.

  for (i=0; i!=m; i++) w[i]=(ARTYPE)0;

  for (i=0; i!=n; i++) {
    t = v[i];
    for (j=pcol[i]; j!=pcol[i+1]; j++) {
      w[irow[j]] += t*a[j];
    }
  }

} // MultMv.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::MultMtv(ARTYPE* v, ARTYPE* w)
{

  int    i,j;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluNonSymMatrix::MultMtv");
  }

  // Determining w = M'.v.

  for (i=0; i!=n; i++) {
    t = (ARTYPE)0;
    for (j=pcol[i]; j!=pcol[i+1]; j++) {
      t += v[irow[j]]*a[j];
    }
    w[i] = t;
  }

} // MultMtv.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::MultMtMv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE* t = new ARTYPE[m];

  MultMv(v,t);
  MultMtv(t,w);

  delete[] t;

} // MultMtMv.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::MultMMtv(ARTYPE* v, ARTYPE* w)
{

  ARTYPE* t = new ARTYPE[n];

  MultMtv(v,t);
  MultMv(t,w);

  delete[] t;

} // MultMMtv.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::Mult0MMt0v(ARTYPE* v, ARTYPE* w)
{

  MultMv(&v[m],w);
  MultMtv(v,&w[m]);

} // Mult0MMt0v.


template<class ARTYPE, class ARFLOAT>
void ARluNonSymMatrix<ARTYPE, ARFLOAT>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARluNonSymMatrix::MultInvv");
  }

  // Solving A.w = v (or AsI.w = v).

  int         info;
  SuperMatrix B;

  if (&v != &w) copy(n, v, 1, w, 1);
  Create_Dense_Matrix(&B, n, 1, w, n, DN, GE);
  gstrs("N", &L, &U, permr, permc, &B, &info);
  Destroy_SuperMatrix_Store(&B); // delete B.Store;

} // MultInvv.


template<class ARTYPE, class ARFLOAT>
inline void ARluNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
             double thresholdp, int orderp, bool check)
{

  m         = np;
  n         = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[n]   = nnz;
  threshold = thresholdp;
  order     = orderp;

  // Checking data.

  if ((check)&&(!DataOK())) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARluSymMatrix::DefineMatrix");
  }

  // Creating SuperMatrix A.

  Create_CompCol_Matrix(&A, n, n, nnz, a, irow, pcol, NC, GE);

  // Reserving memory for vectors used in matrix decomposition.

  permc = new int[n];
  permr = new int[n];

  defined = true;

} // DefineMatrix (square).


template<class ARTYPE, class ARFLOAT>
inline void ARluNonSymMatrix<ARTYPE, ARFLOAT>::
DefineMatrix(int mp, int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp)
{

  m       = mp;
  n       = np;
  nnz     = nnzp;
  a       = ap;
  irow    = irowp;
  pcol    = pcolp;
  pcol[n] = nnz;
  defined = true;
  permc   = NULL;
  permr   = NULL;

} // DefineMatrix (rectangular).


template<class ARTYPE, class ARFLOAT>
inline ARluNonSymMatrix<ARTYPE, ARFLOAT>::ARluNonSymMatrix(): ARMatrix<ARTYPE>()
{ 

  factored = false;  
  permc    = NULL;
  permr    = NULL;

} // Short constructor.


template<class ARTYPE, class ARFLOAT>
inline ARluNonSymMatrix<ARTYPE, ARFLOAT>::
ARluNonSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                 int* pcolp, double thresholdp,
                 int orderp, bool check)                : ARMatrix<ARTYPE>(np)
{

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, thresholdp, orderp, check);

} // Long constructor (square matrix).


template<class ARTYPE, class ARFLOAT>
inline ARluNonSymMatrix<ARTYPE, ARFLOAT>::
ARluNonSymMatrix(int mp, int np, int nnzp, ARTYPE* ap,
                 int* irowp, int* pcolp)            : ARMatrix<ARTYPE>(mp, np)
{

  factored = false;
  DefineMatrix(mp, np, nnzp, ap, irowp, pcolp);

} // Long constructor (retangular matrix).


template<class ARTYPE, class ARFLOAT>
ARluNonSymMatrix<ARTYPE, ARFLOAT>::
ARluNonSymMatrix(char* file, double thresholdp, int orderp, bool check)
{

  factored = false;

  try {
    mat.Define(file);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARluNonSymMatrix");
  }

  if (mat.NCols()==mat.NRows()) {
    DefineMatrix(mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), thresholdp, orderp, check);
  }
  else {                             
    DefineMatrix(mat.NRows(), mat.NCols(), mat.NonZeros(),
                 (ARTYPE*)mat.Entries(), mat.RowInd(), mat.ColPtr());
  }

} // Long constructor (Harwell-Boeing file).


template<class ARTYPE, class ARFLOAT>
ARluNonSymMatrix<ARTYPE, ARFLOAT>& ARluNonSymMatrix<ARTYPE, ARFLOAT>::
operator=(const ARluNonSymMatrix<ARTYPE, ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLNSMAT_H