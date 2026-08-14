using ComplexNumber = System.Numerics.Complex;

using MathNet.Numerics.LinearAlgebra;
using DoubleVector = MathNet.Numerics.LinearAlgebra.Vector<double>;
using DoubleMatrix = MathNet.Numerics.LinearAlgebra.Matrix<double>;

using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

using ComplexLSQSolverResult = LSQSolver.Complex.ComplexLSQSolverResult;

namespace LSQSolver.MathNet
{
    public static class MatrixExtensions
    {
        const double EPS = 2.22044604925032e-16;

        /// <summary>
        /// Solves the least-squares problem Ax ≈ b.
        /// </summary>
        /// <param name="A">The coefficient matrix.</param>
        /// <param name="rhs">The right-hand side vector.</param>
        /// <param name="store_intermediates">If true, stores QR intermediates in the solver result.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrix and vector for NaN and Infinity.</param>
        /// <returns>The available least-squares solution.</returns>
        public static DoubleVector SolveByLSQSolver(this DoubleMatrix A, DoubleVector rhs, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            return SolveByLSQSolver(A, rhs, out _, store_intermediates, rank_tolerance, check_finite);
        }

        /// <summary>
        /// Solves the least-squares problem Ax ≈ b and returns both the solution and the detailed LSQSolver result.
        /// </summary>
        /// <remarks>
        /// If LSQSolver produces a solution with a non-success status, such as
        /// <see cref="LSQSolverStatus.CholeskyFailed"/>, the available solution is returned.
        /// The caller must inspect <see cref="LSQSolverResult.Status"/> before using it.
        /// </remarks>
        /// <param name="A">The coefficient matrix.</param>
        /// <param name="rhs">The right-hand side vector.</param>
        /// <param name="result">The detailed result produced by LSQSolver.</param>
        /// <param name="store_intermediates">If true, stores QR intermediates in <paramref name="result"/>.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrix and vector for NaN and Infinity.</param>
        /// <returns>The available least-squares solution.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="rhs"/> is null.</exception>
        /// <exception cref="ArgumentException">The inputs are empty, incompatible, or contain invalid values.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank_tolerance"/> is negative or non-finite.</exception>
        /// <exception cref="InvalidOperationException">LSQSolver does not produce a solution.</exception>
        public static DoubleVector SolveByLSQSolver(this DoubleMatrix A, DoubleVector rhs, out LSQSolverResult result, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            ArgumentNullException.ThrowIfNull(A);
            ArgumentNullException.ThrowIfNull(rhs);
            ValidateArguments(A.RowCount, A.ColumnCount, rhs.Count, 1, rank_tolerance);

            result = LSQSolver.Solve(
                A.ToColumnMajorArray(), A.RowCount, A.ColumnCount, rhs.ToArray(),
                overwrite: true,
                store_intermediates: store_intermediates,
                rank_tolerance: rank_tolerance,
                check_finite: check_finite);

            ThrowIfNoSolution(result.Status, result.Solution.Length);
            return DoubleVector.Build.DenseOfArray(result.Solution);
        }

        /// <summary>
        /// Solves the multiple-right-hand-side least-squares problem AX ≈ B.
        /// </summary>
        /// <param name="A">The coefficient matrix.</param>
        /// <param name="rhs">The right-hand sides stored as matrix columns.</param>
        /// <param name="store_intermediates">If true, stores QR intermediates in the solver result.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrices for NaN and Infinity.</param>
        /// <returns>The available least-squares solution matrix.</returns>
        public static DoubleMatrix SolveByLSQSolver(this DoubleMatrix A, DoubleMatrix rhs, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            return SolveByLSQSolver(A, rhs, out _, store_intermediates, rank_tolerance, check_finite);
        }

        /// <summary>
        /// Solves the multiple-right-hand-side least-squares problem AX ≈ B and returns the detailed LSQSolver result.
        /// </summary>
        /// <remarks>
        /// Each column of <paramref name="rhs"/> is treated as one right-hand side.
        /// If LSQSolver produces an available solution with a non-success status, that solution is returned.
        /// </remarks>
        /// <param name="A">The coefficient matrix.</param>
        /// <param name="rhs">The right-hand sides stored as matrix columns.</param>
        /// <param name="result">The detailed result produced by LSQSolver.</param>
        /// <param name="store_intermediates">If true, stores QR intermediates in <paramref name="result"/>.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrices for NaN and Infinity.</param>
        /// <returns>The available least-squares solution matrix.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="rhs"/> is null.</exception>
        /// <exception cref="ArgumentException">The inputs are empty, incompatible, or contain invalid values.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank_tolerance"/> is negative or non-finite.</exception>
        /// <exception cref="InvalidOperationException">LSQSolver does not produce a solution.</exception>
        public static DoubleMatrix SolveByLSQSolver(this DoubleMatrix A, DoubleMatrix rhs, out LSQSolverResult result, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            ArgumentNullException.ThrowIfNull(A);
            ArgumentNullException.ThrowIfNull(rhs);
            ValidateArguments(A.RowCount, A.ColumnCount, rhs.RowCount, rhs.ColumnCount, rank_tolerance);

            result = LSQSolver.Solve(
                A.ToColumnMajorArray(), A.RowCount, A.ColumnCount, rhs.ToColumnMajorArray(), rhs.ColumnCount,
                overwrite: true,
                store_intermediates: store_intermediates,
                rank_tolerance: rank_tolerance,
                check_finite: check_finite);

            ThrowIfNoSolution(result.Status, result.Solution.Length);
            return DoubleMatrix.Build.DenseOfColumnMajor(result.Cols, result.RHSCount, result.Solution);
        }

        /// <summary>
        /// Solves the complex least-squares problem Ax ≈ b.
        /// </summary>
        /// <param name="A">The complex coefficient matrix.</param>
        /// <param name="rhs">The complex right-hand side vector.</param>
        /// <param name="store_intermediates">If true, stores intermediates for the realified problem.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrix and vector for NaN and Infinity.</param>
        /// <returns>The available complex least-squares solution.</returns>
        public static ComplexVector SolveByLSQSolver(this ComplexMatrix A, ComplexVector rhs, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            return SolveByLSQSolver(A, rhs, out _, store_intermediates, rank_tolerance, check_finite);
        }

        /// <summary>
        /// Solves the complex least-squares problem Ax ≈ b and returns the detailed LSQSolver.Complex result.
        /// </summary>
        /// <remarks>
        /// The complex problem is realified and solved by the underlying LSQSolver.
        /// If an available solution is produced with a non-success status, such as
        /// <see cref="LSQSolverStatus.CholeskyFailed"/>, that solution is returned.
        /// Intermediate data in <see cref="ComplexLSQSolverResult.KernelResult"/> refers to the realified problem.
        /// </remarks>
        /// <param name="A">The complex coefficient matrix.</param>
        /// <param name="rhs">The complex right-hand side vector.</param>
        /// <param name="result">The detailed result produced by LSQSolver.Complex.</param>
        /// <param name="store_intermediates">If true, stores intermediates for the realified problem.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrix and vector for NaN and Infinity.</param>
        /// <returns>The available complex least-squares solution.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="rhs"/> is null.</exception>
        /// <exception cref="ArgumentException">The inputs are empty, incompatible, or contain invalid values.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank_tolerance"/> is negative or non-finite.</exception>
        /// <exception cref="InvalidOperationException">LSQSolver.Complex does not produce a solution.</exception>
        public static ComplexVector SolveByLSQSolver(this ComplexMatrix A, ComplexVector rhs, out ComplexLSQSolverResult result, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            ArgumentNullException.ThrowIfNull(A);
            ArgumentNullException.ThrowIfNull(rhs);
            ValidateArguments(A.RowCount, A.ColumnCount, rhs.Count, 1, rank_tolerance);

            result = Complex.Solver.Solve(
                A.ToColumnMajorArray(), A.RowCount, A.ColumnCount, rhs.ToArray(),
                store_intermediates: store_intermediates,
                rank_tolerance: rank_tolerance,
                check_finite: check_finite);

            ThrowIfNoSolution(result.Status, result.Solution.Length);
            return ComplexVector.Build.DenseOfArray(result.Solution);
        }

        /// <summary>
        /// Solves the complex multiple-right-hand-side least-squares problem AX ≈ B.
        /// </summary>
        /// <param name="A">The complex coefficient matrix.</param>
        /// <param name="rhs">The complex right-hand sides stored as matrix columns.</param>
        /// <param name="store_intermediates">If true, stores intermediates for the realified problem.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrices for NaN and Infinity.</param>
        /// <returns>The available complex least-squares solution matrix.</returns>
        public static ComplexMatrix SolveByLSQSolver(this ComplexMatrix A, ComplexMatrix rhs, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            return SolveByLSQSolver(A, rhs, out _, store_intermediates, rank_tolerance, check_finite);
        }

        /// <summary>
        /// Solves the complex multiple-right-hand-side least-squares problem AX ≈ B and returns the detailed result.
        /// </summary>
        /// <remarks>
        /// Each column of <paramref name="rhs"/> is treated as one right-hand side.
        /// The complex problem is realified before being passed to the underlying LSQSolver.
        /// Intermediate data in <see cref="ComplexLSQSolverResult.KernelResult"/> refers to the realified problem.
        /// </remarks>
        /// <param name="A">The complex coefficient matrix.</param>
        /// <param name="rhs">The complex right-hand sides stored as matrix columns.</param>
        /// <param name="result">The detailed result produced by LSQSolver.Complex.</param>
        /// <param name="store_intermediates">If true, stores intermediates for the realified problem.</param>
        /// <param name="rank_tolerance">The relative tolerance used for numerical rank detection.</param>
        /// <param name="check_finite">If true, checks the input matrices for NaN and Infinity.</param>
        /// <returns>The available complex least-squares solution matrix.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="rhs"/> is null.</exception>
        /// <exception cref="ArgumentException">The inputs are empty, incompatible, or contain invalid values.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="rank_tolerance"/> is negative or non-finite.</exception>
        /// <exception cref="InvalidOperationException">LSQSolver.Complex does not produce a solution.</exception>
        public static ComplexMatrix SolveByLSQSolver(this ComplexMatrix A, ComplexMatrix rhs, out ComplexLSQSolverResult result, bool store_intermediates = false, double rank_tolerance = EPS, bool check_finite = true)
        {
            ArgumentNullException.ThrowIfNull(A);
            ArgumentNullException.ThrowIfNull(rhs);
            ValidateArguments(A.RowCount, A.ColumnCount, rhs.RowCount, rhs.ColumnCount, rank_tolerance);

            result = Complex.Solver.Solve(
                A.ToColumnMajorArray(), A.RowCount, A.ColumnCount, rhs.ToColumnMajorArray(), rhs.ColumnCount,
                store_intermediates: store_intermediates,
                rank_tolerance: rank_tolerance,
                check_finite: check_finite);

            ThrowIfNoSolution(result.Status, result.Solution.Length);
            return ComplexMatrix.Build.DenseOfColumnMajor(result.Cols, result.RHSCount, result.Solution);
        }

        /// <summary>
        /// Validates matrix dimensions, right-hand-side dimensions, and rank tolerance.
        /// </summary>
        private static void ValidateArguments(int rows, int cols, int rhsRows, int rhsCount, double rank_tolerance)
        {
            if (rows <= 0 || cols <= 0)
                throw new ArgumentException("The coefficient matrix must not be empty.");

            if (rhsRows != rows)
                throw new ArgumentException("The number of right-hand-side rows must equal the number of coefficient matrix rows.");

            if (rhsCount <= 0)
                throw new ArgumentException("At least one right-hand side is required.");

            if (!double.IsFinite(rank_tolerance) || rank_tolerance < 0.0)
                throw new ArgumentOutOfRangeException(nameof(rank_tolerance), "The rank tolerance must be finite and non-negative.");
        }

        /// <summary>
        /// Throws when the solver does not provide an available solution.
        /// </summary>
        private static void ThrowIfNoSolution(LSQSolverStatus status, int solutionLength)
        {
            if (solutionLength > 0) return;

            if (status is LSQSolverStatus.NullMatrix
                or LSQSolverStatus.EmptyMatrix
                or LSQSolverStatus.NullVector
                or LSQSolverStatus.DimensionMismatch
                or LSQSolverStatus.InvalidVector
                or LSQSolverStatus.InvalidMatrix
                or LSQSolverStatus.InvalidMatrixStorage)
                throw new ArgumentException($"LSQSolver rejected the input. Status: {status}.");

            throw new InvalidOperationException($"LSQSolver did not produce a solution. Status: {status}.");
        }
    }
}