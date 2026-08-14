# Changelog

All notable changes to `LSQSolver.MathNet` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-08-14

### Added

- Initial release of `LSQSolver.MathNet`.
- Added `SolveByLSQSolver(...)` extension methods for MathNet dense matrices.
- Added real-valued least-squares support for:
  - `Matrix<double>` with a `Vector<double>` right-hand side.
  - `Matrix<double>` with a `Matrix<double>` right-hand side.
- Added complex-valued least-squares support for:
  - `Matrix<System.Numerics.Complex>` with a `Vector<System.Numerics.Complex>` right-hand side.
  - `Matrix<System.Numerics.Complex>` with a `Matrix<System.Numerics.Complex>` right-hand side.
- Added support for overdetermined, square, underdetermined, and numerically rank-deficient systems through the same extension-method API.
- Added minimum-2-norm least-squares solutions when the solution is not unique.
- Added multiple-right-hand-side solves, with each column of the right-hand-side matrix treated as one right-hand side.
- Added overloads that return only the available MathNet solution vector or matrix.
- Added overloads with an `out` parameter that expose the detailed solver result:
  - `LSQSolverResult` for real problems.
  - `ComplexLSQSolverResult` for complex problems.
- Added optional control of:
  - QR-intermediate storage through `store_intermediates`.
  - Numerical-rank detection through `rank_tolerance`.
  - Input validation for `NaN` and infinity through `check_finite`.
- Added argument validation for empty inputs, incompatible dimensions, invalid right-hand-side counts, and invalid rank tolerances.
- Added exception handling that distinguishes invalid arguments from cases in which the solver produces no available solution.
- Added support for returning an available solution together with a non-success solver status, allowing callers to inspect the status before deciding whether to use the solution.

### Implementation Notes

- MathNet matrices and vectors are copied to column-major work arrays before calling LSQSolver.
- The LSQSolver kernel is invoked in overwrite mode on those private work arrays; the original MathNet inputs are not modified.
- Complex systems are converted to equivalent real least-squares systems and solved by `LSQSolver.Complex` using the real LSQSolver kernel.
- For complex problems, intermediate data in `ComplexLSQSolverResult.KernelResult` refers to the realified problem.

### Documentation

- Added English and Japanese README files.
- Added usage examples for:
  - Basic real-valued solves.
  - Diagnostic result handling.
  - Underdetermined and rank-deficient problems.
  - Multiple right-hand sides.
  - Complex-valued problems.
  - Optional solver parameters and error handling.
- Documented the intended distinction between MathNet's standard `Solve()` and `SolveByLSQSolver()`.
- Documented `LSQSolver.MathNet` as one possible external solution for underdetermined least-squares use cases discussed in MathNet.Numerics issue [#560](https://github.com/mathnet/mathnet-numerics/issues/560).
- Added references to related MathNet.Numerics discussions [#490](https://github.com/mathnet/mathnet-numerics/issues/490) and [#580](https://github.com/mathnet/mathnet-numerics/issues/580), while clarifying that this package does not directly fix those issues.
