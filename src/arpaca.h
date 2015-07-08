#ifndef ARPACA_ARPACA_H
#define ARPACA_ARPACA_H

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <eigen3/Eigen/Core>

// ARPACK interface
extern "C" {

// Lanczos method (double precision)
void dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
             double* resid, int* ncv, double* v, int* ldv, int* iparam,
             int* ipntr, double* workd, double* workl, int* lworkl, int* info);

// Lanczos method (single precision)
void ssaupd_(int* ido, char* bmat, int* n, char* which, int* nev, float* tol,
             float* resid, int* ncv, float* v, int* ldv, int* iparam,
             int* ipntr, float* workd, float* workl, int* lworkl, int* info);

// Eigenvalue problem solver for result of dsaupd_
void dseupd_(int* rvec, char* howmny, int* select, double* d, double* z,
             int* ldz, double* sigma, char* bmat, int* n, char* which, int* nev,
             double* tol, double* resid, int* ncv, double* v, int* ldv,
             int* iparam, int* ipntr, double* workd, double* workl, int* lworkl,
             int* info);

// Eigenvalue problem solver for result of ssaupd_
void sseupd_(int* rvec, char* howmny, int* select, float* d, float* z, int* ldz,
             float* sigma, char* bmat, int* n, char* which, int* nev,
             float* tol, float* resid, int* ncv, float* v, int* ldv,
             int* iparam, int* ipntr, float* workd, float* workl, int* lworkl,
             int* info);

}  // extern "C"

namespace arpaca {

// Type of eigenvalue to solve.
// Arpaca::SymmetricEigenSolver computes part of eigenvalues and eigenvectors.
// EigenvalueType indicates which part is to be solved.
enum EigenvalueType {
    ALGEBRAIC_LARGEST,  // Largest as real numbers
    ALGEBRAIC_SMALLEST,  // Smallest as real numbers
    ALGEBRAIC_BOTH_END,  // Both end as real numbers
    MAGNITUDE_LARGEST,  // Largest as absolute values
    MAGNITUDE_SMALLEST,  // Smallest as absolute values
};

namespace detail {

inline const char* GetNameOfEigenvalueType(EigenvalueType type)
{
    switch (type) {
    case ALGEBRAIC_LARGEST: return "LA";
    case ALGEBRAIC_SMALLEST: return "SA";
    case ALGEBRAIC_BOTH_END: return "BE";
    case MAGNITUDE_LARGEST: return "LM";
    case MAGNITUDE_SMALLEST: return "SM";
    }
    throw std::invalid_argument("Invalid eigenvalue type");
}

inline void GetNameOfEigenvalueType(EigenvalueType type, char* to)
{
    const char* from = GetNameOfEigenvalueType(type);
    std::copy(from, from + 3, to);
}

template<typename Scalar>
void saupd(int* ido, char* bmat, int* n, char* which, int* nev, Scalar* tol,
           Scalar* resid, int* ncv, Scalar* v, int* ldv, int* iparam,
           int* ipntr, Scalar* workd, Scalar* workl, int* lworkl, int* info)
{
    dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr,
            workd, workl, lworkl, info);
}

template<>
inline void saupd<float>(int* ido, char* bmat, int* n, char* which, int* nev,
                         float* tol, float* resid, int* ncv, float* v, int* ldv,
                         int* iparam, int* ipntr, float* workd, float* workl,
                         int* lworkl, int* info)
{
    ssaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr,
            workd, workl, lworkl, info);
}

template<typename Scalar>
void seupd(int* rvec, char* howmny, int* select, Scalar* d, Scalar* z, int* ldz,
           Scalar* sigma, char* bmat, int* n, char* which, int* nev,
           Scalar* tol, Scalar* resid, int* ncv, Scalar* v, int* ldv,
           int* iparam, int* ipntr, Scalar* workd, Scalar* workl, int* lworkl,
           int* info)
{
    dseupd_(rvec, howmny, select, d, z, ldz, sigma, bmat, n, which, nev, tol,
            resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

template<>
inline void seupd<float>(int* rvec, char* howmny, int* select, float* d,
                         float* z, int* ldz, float* sigma, char* bmat, int* n,
                         char* which, int* nev, float* tol, float* resid,
                         int* ncv, float* v, int* ldv, int* iparam, int* ipntr,
                         float* workd, float* workl, int* lworkl, int* info)
{
    sseupd_(rvec, howmny, select, d, z, ldz, sigma, bmat, n, which, nev, tol,
            resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

}  // namespace detail

// Eigenvalue problem solver for symmetric matrix (especially sparse one).
// It partly computes eigenvalues and eigenvectors of a symmetric matrix A using
// ARPACK, which implements Implicitly Restarted Arnoldi (IRA) method.
template<typename Scalar = double>
class SymmetricEigenSolver {
public:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    SymmetricEigenSolver()
        : num_lanczos_vectors_(100),
          max_iterations_(100),
          tolerance_(-1)
    {
        SetEigenvalueType(MAGNITUDE_LARGEST);
    }

    // Sets # of Lanczos vectors.
    // From two to four times larger than # of objective eigenvalues is
    // recommended.
    //
    // NOTE: IRA approximates subspace containing the objective eigenvectors by
    // constructing |num| orthogonal basis vectors in each iteration, which is
    // improved for each iteration. Consequently, too small |num| causes slow
    // convergence, while too large |num| requires large size of memory and even
    // make each iteration too slow.
    void SetNumLanczosVectors(int num)
    {
        num_lanczos_vectors_ = num;
    }

    // Sets max # of Lanczos iterations.
    void SetMaxIterations(int num)
    {
        if (num <= 0)
            throw std::invalid_argument("Given # of iterations is not positive");
        max_iterations_ = num;
    }

    // Sets precision tolerance for computation.
    // Negative value means machine precision.
    void SetTolerance(Scalar tolerance)
    {
        tolerance_ = tolerance;
    }

    // Sets the type of eigenvalues to solve.
    // See EigenvalueType for detail.
    void SetEigenvalueType(EigenvalueType type)
    {
        detail::GetNameOfEigenvalueType(type, eigenvalue_type_name_);
    }

    // Solves the eigenvalue problem of a real symmetrix matrix A.
    // The matrix A is not set explicitly. Solve() only needs the operator
    // y = A * x.
    //
    // dimension:       Dimension of space (rows or cols of the matrix A).
    // num_eigenvalues: # of eigenvalues to compute.
    // op:              Functor that computes operation y = A * x. Signature is
    //                  op(x, y), where x and y are slices of Eigen::VectorX{d,f},
    //                  i.e. return values of Eigen::VectorX{d,f}::middleRows().
    template<typename OperatorATimes>
    void Solve(int dimension, int num_eigenvalues, OperatorATimes op)
    {
        if (dimension <= 0)
            throw std::invalid_argument("Given dimension is not positive");
        if (num_eigenvalues <= 0)
            throw std::invalid_argument("Given num_eigenvalues is not positive");

        // Eigenvalues (Lanczos)
        int info = 0;

        num_eigenvalues = std::max(std::min(num_eigenvalues, dimension), 0);

        int num_lanczos_vectors = num_lanczos_vectors_;
        if (num_lanczos_vectors <= 0)
            num_lanczos_vectors = num_eigenvalues * 3;
        else if (num_lanczos_vectors <= num_eigenvalues)
            num_lanczos_vectors = num_eigenvalues + 1;

        num_lanczos_vectors = std::min(num_lanczos_vectors, dimension);

        char bmat[2] = { 'I' };
        int leading_dimension = dimension;

        int iparam[11] = {}, ipntr[11] = {};
        iparam[0] = 1;  // Use effective shift
        iparam[2] = max_iterations_;
        iparam[3] = 1;  // block size
        iparam[6] = 1;  // Solve standard eigenvalue problem

        // Calculation workspace
        int lworkl = (num_lanczos_vectors + 8) * num_lanczos_vectors;
        Vector workd(3 * dimension), workl(lworkl);
        std::vector<Scalar> residue(dimension);
        Matrix lanczos_vectors(dimension, num_lanczos_vectors);

        // Lanczos iterations
        for (int mode = 0; mode != 99; ) {
            detail::saupd(&mode, bmat, &dimension, eigenvalue_type_name_,
                          &num_eigenvalues, &tolerance_, &residue[0],
                    &num_lanczos_vectors, lanczos_vectors.data(),
                    &leading_dimension, iparam, ipntr, workd.data(),
                    workl.data(), &lworkl, &info);

            switch (mode) {
            case -1:
                op(workd.middleRows(ipntr[0] - 1, dimension),
                        workd.middleRows(ipntr[1] - 1, dimension));
                break;
            case 1:
                op(workd.middleRows(ipntr[2] - 1, dimension),
                        workd.middleRows(ipntr[1] - 1, dimension));
                workd.middleRows(ipntr[2] - 1, dimension) =
                        workd.middleRows(ipntr[0] - 1, dimension);
                break;
            case 2:
                workd.middleRows(ipntr[1] - 1, dimension) =
                        workd.middleRows(ipntr[0] - 1, dimension);
                break;
            }
        }

        info_ = info;
        num_actual_iterations_ = iparam[2];
        num_converged_eigenvalues_ = iparam[4];


        // Eigenvectors
        int rvec = 1;  // Compute eigenvectors
        char howmany[2] = { 'A' };
        std::vector<int> select(num_lanczos_vectors);
        Scalar sigma;

        eigenvectors_.resize(dimension, num_eigenvalues);
        eigenvalues_.resize(num_eigenvalues);

        detail::seupd(&rvec, howmany, &select[0], eigenvalues_.data(),
                eigenvectors_.data(), &dimension, &sigma, bmat, &dimension,
                eigenvalue_type_name_, &num_eigenvalues, &tolerance_,
                &residue[0], &num_lanczos_vectors, lanczos_vectors.data(),
                &leading_dimension, iparam, ipntr, workd.data(), workl.data(),
                &lworkl, &info);
    }

    // Gets the computed eigenvectors.
    // Returned matrix contains each eigenvector in each column.
    const Matrix& eigenvectors() const
    {
        return eigenvectors_;
    }

    // Gets the computed eigenvalues.
    const Vector& eigenvalues() const
    {
        return eigenvalues_;
    }

    Matrix MoveEigenvectors()
    {
        Matrix eigenvectors;
        eigenvectors.swap(eigenvectors_);
        return eigenvectors;
    }

    Vector MoveEigenvalues()
    {
        Vector eigenvalues;
        eigenvalues.swap(eigenvalues_);
        return eigenvalues;
    }

    int num_actual_iterations() const
    {
        return num_actual_iterations_;
    }

    int num_converged_eigenvalues() const
    {
        return num_converged_eigenvalues_;
    }

    // Gets information message based on error code of ARPACK.
    const char* GetInfo() const
    {
        switch (info_) {
        case 0: return "Normal exit";
        case 1: return "Not converged";
        case 3: return "num_lanczos_vectors may be too small";
        case -8: return "Error on tridiagonal eigenvalue calculation";
        case -9999: return "Could not build an Arnoldi factorization";
        }
        return "Arpaca internal error";
    }

private:
    Matrix eigenvectors_;
    Vector eigenvalues_;

    int num_lanczos_vectors_;
    int max_iterations_;
    Scalar tolerance_;

    int info_;
    int num_actual_iterations_;
    int num_converged_eigenvalues_;

    char eigenvalue_type_name_[3];
};

// Default operator for the case the matrix A is explicitly given.
template<typename MatrixType>
class DefaultOperator {
public:
    explicit DefaultOperator(MatrixType& m)
        : m_(m)
    {}

    template<typename X, typename Y>
    void operator()(X x, Y y) const
    {
        y = m_ * x;
    }

private:
    MatrixType& m_;
};

template<typename Matrix>
inline DefaultOperator<Matrix> MakeDefaultOperator(Matrix& A)
{
    return DefaultOperator<Matrix>(A);
}

// Solves the eigenvalue problem of a real symmetrix matrix A, where the matrix
// A is explicitly given.
template<typename Matrix>
SymmetricEigenSolver<typename Matrix::Scalar> Solve(
        const Matrix& A,
        int num_eigenvalues,
        EigenvalueType type = MAGNITUDE_LARGEST,
        int num_lanczos_vectors = 0,
        int max_iterations = 10000,
        typename Matrix::Scalar tolerance = -1)
{
    SymmetricEigenSolver<typename Matrix::Scalar> solver;
    solver.SetEigenvalueType(type);
    solver.SetMaxIterations(max_iterations);
    solver.SetNumLanczosVectors(num_lanczos_vectors);
    solver.SetTolerance(tolerance);

    solver.Solve(A.rows(), num_eigenvalues, MakeDefaultOperator(A));

    return solver;
}

}  // namespace arpaca

#endif  // ARPACA_ARPACA_HPP_
