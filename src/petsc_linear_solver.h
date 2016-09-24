#ifndef PETSC_LINEAR_SOLVER_H
#define PETSC_LINEAR_SOLVER_H

#include <petscksp.h>
#include <zjucad/matrix/matrix.h>

#ifndef PETSC_VERSION_GT
#define PETSC_VERSION_(MAJOR,MINOR,SUBMINOR)    \
  ((PETSC_VERSION_MAJOR == (MAJOR)) &&          \
   (PETSC_VERSION_MINOR == (MINOR)) &&          \
   (PETSC_VERSION_SUBMINOR == (SUBMINOR)) &&    \
   (PETSC_VERSION_RELEASE  == 1))

#define PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR)          \
  (PETSC_VERSION_RELEASE == 1 &&                        \
   (PETSC_VERSION_MAJOR < (MAJOR) ||                    \
    (PETSC_VERSION_MAJOR == (MAJOR) &&                  \
     (PETSC_VERSION_MINOR < (MINOR) ||                  \
      (PETSC_VERSION_MINOR == (MINOR) &&                \
       (PETSC_VERSION_SUBMINOR < (SUBMINOR)))))))

#define PETSC_VERSION_LE(MAJOR,MINOR,SUBMINOR)  \
  (PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR) ||    \
   PETSC_VERSION_(MAJOR,MINOR,SUBMINOR))

#define PETSC_VERSION_GT(MAJOR,MINOR,SUBMINOR)  \
  (!PETSC_VERSION_LE(MAJOR,MINOR,SUBMINOR))
#endif


#if PETSC_VERSION_GT(3, 1, 0)
#define vec_destroy(v) VecDestroy(v)
#define mat_destroy(m) MatDestroy(m)
#define KSP_destroy(k) KSPDestroy(k)
#else
#define vec_destroy(v) VecDestroy(*v)
#define mat_destroy(m) MatDestroy(*m)
#define KSP_destroy(k) KSPDestroy(*k)
#endif

#if PETSC_VERSION_GT(3, 2, 0)
#define vec_create_seq_with_array(comm, bs, n, a, v)    \
  VecCreateSeqWithArray(comm, bs, n, a, v)
#else
#define vec_create_seq_with_array(comm, bs, n, a, v)    \
  VecCreateSeqWithArray(comm, n, a, v)
#endif

#if PETSC_VERSION_LT(3, 5, 0)
#  define KSPSetOperators_SAME_PC(ksp, A, B) KSPSetOperators(ksp, A, B, SAME_PRECONDITIONER)
#else
#  define KSPSetOperators_SAME_PC(ksp, A, B) \
  KSPSetOperators(ksp, A, B); \
  KSPSetReusePreconditioner(ksp, PETSC_TRUE)
#endif

#if PETSC_VERSION_GT(3, 4, 0)
  #include <petsc/private/matimpl.h>
#else
  #if PETSC_VERSION_GT(3, 2, 0)
    #include <petsc-private/matimpl.h>
  #endif
  #include <private/matimpl.h>
#endif

#include <hjlib/sparse/sparse.h>
#include <zjucad/matrix/matrix.h>

class PETsc_imp
{
 public:
  PETsc_imp() {
    std::cout << "call PETsc_imp" << std::endl;
    if(!PetscInitializeCalled){
      int argc = 0;
      PetscInitializeNoArguments();
    }
  }
  virtual ~PETsc_imp() {
    if(!PetscFinalizeCalled){
      PetscFinalize();
    }
  }
};

class PETsc_CG_imp
{
 public:
  PETsc_CG_imp(const double * val, const int32_t * idx,  const int32_t * ptr,
               const size_t nnz, const size_t row, const size_t col, const char *pc_str) {
    std::cout << "call PETsc_CG_imp" << std::endl;
    comm = MPI_COMM_SELF;
    Dim = row;
    if(sizeof(int32_t) == sizeof(int)) {
      MatCreateSeqAIJWithArrays(
          comm,
          row, col,
          const_cast<int *>(ptr), const_cast<int *>(idx), const_cast<double *>(val), &A_);
    }
    else {
      zjucad::matrix::itr_matrix<const int32_t*> ptr_m(col, 1, ptr);
      zjucad::matrix::itr_matrix<const int32_t*> idx_m(nnz, 1, idx);
      ptr_ = ptr_m;
      idx_ = idx_m;
      MatCreateSeqAIJWithArrays(
          comm,
          row, col, &ptr_[0], &idx_[0], const_cast<double *>(val), &A_);
    }
    KSPCreate(comm,&solver);
    KSPSetOperators_SAME_PC(solver, A_, A_);
    PC pc;
    KSPGetPC(solver,&pc);
    PCSetType(pc, pc_str);
    KSPSetType(solver, KSPCG);
    KSPCGSetType(solver, KSP_CG_SYMMETRIC);
    KSPSetInitialGuessNonzero(solver,PETSC_TRUE);
    KSPSetUp(solver);
  }
  int solve(const double *b, double *x, size_t rhs) {
    KSPSetTolerances(solver,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    vec_create_seq_with_array(comm, 1, Dim, b, &B_);
    vec_create_seq_with_array(comm, 1, Dim, x, &X_);
    KSPSolve(solver, B_, X_);

    vec_destroy(&B_);
    vec_destroy(&X_);
  }
  ~PETsc_CG_imp() {
    mat_destroy(&A_);
    KSP_destroy(&solver);
  }
private:
  Mat A_;
  Vec X_, B_;
  zjucad::matrix::matrix<int> ptr_, idx_;
  KSP solver;
  MPI_Comm comm;
  int Dim;
};
  
#endif
