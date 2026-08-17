# LSQSolver.MathNet

[日本語版](https://github.com/TaigaNakano/LSQSolver.MathNet/blob/master/README_jp.md)

`LSQSolver.MathNet` adds MathNet.Numerics extension methods for solving dense real and complex least-squares problems with [LSQSolver](https://github.com/TaigaNakano/LSQSolver).

It is intended especially for cases that are not naturally covered by the standard `Matrix<T>.Solve(...)` API:

- underdetermined systems (`rows < columns`),
- rank-deficient systems,
- least-squares problems requiring a minimum-2-norm solution, and
- applications that need numerical-rank and residual diagnostics.

For regular square systems, MathNet.Numerics' optimized `Solve(...)` is usually the better choice and may be substantially faster. This package provides consistent least-squares semantics across matrix shapes and numerical ranks; it is not a replacement for MathNet's optimized LU solver.

**Related projects:** [LSQSolver](https://github.com/TaigaNakano/LSQSolver) · [LSQSolver.Complex](https://github.com/TaigaNakano/LSQSolver.Complex) · [LSQSolver.MathNet](https://github.com/TaigaNakano/LSQSolver.MathNet)

## Installation

```bash
dotnet add package LSQSolver.MathNet
```

The package uses MathNet.Numerics matrix and vector types and delegates the numerical solve to LSQSolver.

## Basic usage

### Real matrices

```csharp
using LSQSolver.MathNet;
using MathNet.Numerics.LinearAlgebra;

Matrix<double> A = Matrix<double>.Build.DenseOfArray(new[,]
{
    { 1.0, 0.0, 1.0 },
    { 0.0, 1.0, 1.0 }
});

Vector<double> b = Vector<double>.Build.Dense(new[] { 1.0, 1.0 });

Vector<double> x = A.SolveByLSQSolver(b);
```

Here, `A` is underdetermined. `SolveByLSQSolver` returns a least-squares solution with minimum Euclidean norm when the solve succeeds.

### Complex matrices

```csharp
using LSQSolver.MathNet;
using MathNet.Numerics.LinearAlgebra;
using Complex = System.Numerics.Complex;

Matrix<Complex> A = Matrix<Complex>.Build.DenseOfArray(new[,]
{
    { new Complex(1.0, 1.0), Complex.Zero },
    { Complex.One,              Complex.One }
});

Vector<Complex> b = Vector<Complex>.Build.Dense(new[]
{
    Complex.One,
    new Complex(0.0, 1.0)
});

Vector<Complex> x = A.SolveByLSQSolver(b);
```

The complex adapter converts the problem to an equivalent real least-squares system and reconstructs the complex solution.

## Diagnostics

Use the overload with an `out` parameter when the numerical rank, residual norm, solver status, or optional intermediate data is required.

```csharp
Vector<double> x = A.SolveByLSQSolver(
    b,
    out var result,
    store_intermediates: true);

if (result.Status != LSQSolverStatus.Success)
{
    Console.WriteLine($"Solver status: {result.Status}");
}

Console.WriteLine($"Rank: {result.Rank}");
Console.WriteLine($"Residual norm: {result.ResidualNorm}");
Console.WriteLine(result.ToString(
    omit: false,
    display_row_count: 10,
    display_col_count: 10));
```

For complex problems, the adapter result exposes the underlying real solver result through `KernelResult`.

```csharp
Vector<Complex> x = A.SolveByLSQSolver(b, out var result);
Console.WriteLine(result.KernelResult?.Status);
```

## Options

The extension methods provide the following optional arguments:

| Argument | Description |
| --- | --- |
| `store_intermediates` | Stores QR-related intermediate data in the result when `true`. |
| `rank_tolerance` | Relative tolerance used for numerical-rank detection. |
| `check_finite` | Checks the input for `NaN` and infinity when `true`. |

Input MathNet matrices and vectors are not overwritten. They are converted to the column-major arrays required by LSQSolver.

## Solution semantics

LSQSolver computes a solution of

$$
\min_x \|Ax-b\|_2.
$$

If the minimizer is not unique, the solver selects a minimum-2-norm solution:

$$
\min \{\|x\|_2 : x \in \mbox{argmin}_y \|Ay-b\|_2\}.
$$

The same interpretation is used for overdetermined, underdetermined, and numerically rank-deficient systems.

The numerical method is based on column-pivoted QR factorization, numerical-rank detection, and minimum-norm completion. It does not compute a full SVD.

## Relationship to MathNet.Numerics

MathNet.Numerics provides efficient direct solvers for regular square systems and QR-based least-squares solvers for supported rectangular systems. Its standard `Matrix<T>.Solve(...)` API, however, does not naturally cover every underdetermined or rank-deficient least-squares problem.

Related limitations and use cases have been discussed in MathNet.Numerics issues:

- [#560: Cannot solve linear system if input matrix has less rows than columns](https://github.com/mathnet/mathnet-numerics/issues/560)
- [#490: QRFactor error in native providers](https://github.com/mathnet/mathnet-numerics/issues/490)
- [#580: Matrix Inverse NaN/Infinity/-Infinity](https://github.com/mathnet/mathnet-numerics/issues/580)

`LSQSolver.MathNet` is one possible external solution for the least-squares use cases represented most directly by issue #560. It lets existing MathNet matrices call a solver that supports underdetermined and rank-deficient systems without requiring users to select and combine separate factorization APIs themselves.

Issue #580 concerns matrix inversion rather than least-squares solving. This package does not define an inverse for a singular matrix. It is relevant only when the actual goal is to solve or approximate `Ax = b`, in which case a least-squares or minimum-norm solution may be the appropriate operation instead of forming `A.Inverse()`.

MathNet.Numerics also provides SVD and `PseudoInverse()` as explicit alternatives. This package offers a different algorithm and a `Solve`-style interface specialized for dense least-squares problems.

## Choosing between MathNet `Solve` and LSQSolver

| Problem | Suggested approach |
| --- | --- |
| Regular square system | Prefer MathNet `A.Solve(b)` for its optimized LU factorization. |
| Full-column-rank overdetermined system | Either solver may be appropriate; benchmark the actual workload. |
| Underdetermined system | Use `SolveByLSQSolver` when a minimum-2-norm solution is required. |
| Rank-deficient least-squares system | Use `SolveByLSQSolver` when rank-aware minimum-norm handling is required. |
| Explicit pseudoinverse required |For pseudoinverses, also consider MathNet's `PseudoInverse()`. This package can construct one using `A.SolveByLSQSolver(Matrix<double>.Build.DenseIdentity(A.RowCount))`, but compare its performance and accuracy for your problem. When solving `Ax=b`, explicitly forming the pseudoinverse is slower and requires more memory.
|

## Numerical considerations

- Numerical rank depends on the scale of the matrix and `rank_tolerance`.
- A minimum-norm solution is a mathematical selection rule, not necessarily the appropriate physical prior for an inverse problem.
- Severe scaling or conditioning problems may require normalization, regularization, or an SVD-based method.
- Always inspect the returned status before relying on a diagnostic result.

## Related Projects

| Project | Description |
|---|---|
| [LSQSolver](https://github.com/TaigaNakano/LSQSolver) | The core rank-aware least-squares solver for real-valued dense problems. |
| [LSQSolver.Complex](https://github.com/TaigaNakano/LSQSolver.Complex) | Complex-valued least-squares support built on LSQSolver. |
| [LSQSolver.MathNet](https://github.com/TaigaNakano/LSQSolver.MathNet) | MathNet.Numerics integration for real and complex least-squares problems. |

## License

MIT License
