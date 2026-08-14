# LSQSolver.MathNet

[English version](LSQSolver.MathNet-README.md)

`LSQSolver.MathNet`は、[LSQSolver](https://github.com/TaigaNakano/LSQSolver)をMathNet.Numericsの行列・ベクトルから利用するための拡張メソッドを提供します。実数および複素数の密行列最小二乗問題に対応します。

特に、MathNet.Numericsの標準的な`Matrix<T>.Solve(...)`では自然に扱いにくい、次の問題を対象としています。

- 劣決定問題（`rows < columns`）
- ランク落ち問題
- 最小2ノルム解を必要とする最小二乗問題
- 数値ランクや残差などの診断情報を必要とする問題

正則な正方行列では、通常、MathNet.Numericsの最適化された`Solve(...)`を使用する方が適切であり、大幅に高速な場合があります。このパッケージの目的はMathNetのLUソルバーを置き換えることではなく、行列の形状や数値ランクによらず一貫した意味で最小二乗問題を扱えるようにすることです。

## インストール

```bash
dotnet add package LSQSolver.MathNet
```

このパッケージはMathNet.Numericsの行列・ベクトル型を受け取り、数値計算をLSQSolverへ委譲します。

## 基本的な使用方法

### 実数行列

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

この例の`A`は劣決定行列です。解法が成功した場合、`SolveByLSQSolver`は最小二乗解のうちユークリッドノルムが最小の解を返します。

### 複素行列

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

複素数アダプターは、複素最小二乗問題を等価な実最小二乗問題へ変換し、計算後に複素数解を再構成します。

## 診断情報

数値ランク、残差ノルム、ソルバーの状態、または中間情報が必要な場合は、`out`引数を持つオーバーロードを使用します。

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

複素数問題では、アダプターの結果から`KernelResult`を介して内部の実数ソルバー結果を参照できます。

```csharp
Vector<Complex> x = A.SolveByLSQSolver(b, out var result);
Console.WriteLine(result.KernelResult?.Status);
```

## オプション

拡張メソッドは次のオプション引数を持ちます。

| 引数 | 説明 |
| --- | --- |
| `store_intermediates` | `true`の場合、QR分解に関係する中間情報を結果へ保存します。 |
| `rank_tolerance` | 数値ランク判定に使用する相対許容値です。 |
| `check_finite` | `true`の場合、入力に`NaN`または無限大が含まれていないか検査します。 |

入力されたMathNetの行列とベクトルは上書きされません。LSQSolverが必要とする列優先配列へ変換してから計算します。

## 解の意味

LSQSolverは次の問題の解を計算します。

$$
\min_x \|Ax-b\|_2.
$$

最小化解が一意でない場合、その中から最小2ノルム解を選択します。

$$
\min \{\|x\|_2 : x \in \mbox{argmin}_y \|Ay-b\|_2\}.
$$

過決定、劣決定、および数値的にランク落ちした問題を、同じ意味で取り扱います。

数値解法は列ピボット付きQR分解、数値ランク判定、および最小ノルム補完に基づきます。完全なSVDは計算しません。

## MathNet.Numericsとの関係

MathNet.Numericsは、正則な正方行列に対する効率的な直接法と、対応する長方形行列に対するQRベースの最小二乗解法を提供しています。一方、標準的な`Matrix<T>.Solve(...)` APIだけでは、すべての劣決定問題やランク落ち最小二乗問題を自然に扱えるわけではありません。

関連する制約やユースケースは、MathNet.NumericsのIssueでも議論されています。

- [#560: Cannot solve linear system if input matrix has less rows than columns](https://github.com/mathnet/mathnet-numerics/issues/560)
- [#490: QRFactor error in native providers](https://github.com/mathnet/mathnet-numerics/issues/490)
- [#580: Matrix Inverse NaN/Infinity/-Infinity](https://github.com/mathnet/mathnet-numerics/issues/580)

`LSQSolver.MathNet`は、特にIssue #560に示されている最小二乗問題に対する外部Solutionの一つです。既存のMathNet行列から、劣決定およびランク落ち問題に対応するソルバーを直接呼び出せるため、利用者が複数の分解APIを個別に選択・組み合わせる必要を減らします。

Issue #580は最小二乗問題ではなく逆行列計算に関するものです。このパッケージは、特異行列の逆行列を定義するものではありません。実際の目的が`Ax = b`を解く、または近似的に解くことである場合に限り、`A.Inverse()`を構成する代わりに最小二乗解や最小ノルム解を求める方法が関連します。

MathNet.Numericsには、明示的な代替手段としてSVDおよび`PseudoInverse()`も用意されています。このパッケージは、密な最小二乗問題に特化した別のアルゴリズムと`Solve`形式のインターフェースを提供します。

## MathNetの`Solve`との使い分け

| 問題 | 推奨される方法 |
| --- | --- |
| 正則な正方行列 | 最適化されたLU経路を使用するMathNetの`A.Solve(b)`を優先します。 |
| 列フルランクの過決定問題 | どちらも候補になります。実際の問題で性能と精度を比較してください。 |
| 劣決定問題 | 最小2ノルム解が必要な場合は`SolveByLSQSolver`を使用します。 |
| ランク落ち最小二乗問題 | ランクを考慮した最小ノルム処理が必要な場合は`SolveByLSQSolver`を使用します。 |
| 擬似逆行列そのものが必要 | MathNetの`PseudoInverse()`を検討してください。本パッケージを逆行列演算として使用しないでください。 |

## 数値計算上の注意

- 数値ランクは、行列のスケールと`rank_tolerance`に依存します。
- 最小ノルム解は数学的な解選択規則であり、逆問題に対する物理的な事前情報として常に適切とは限りません。
- 著しいスケーリングまたは悪条件性がある場合、正規化、正則化、あるいはSVDベースの解法が必要になることがあります。
- 診断結果を利用する場合は、必ず返されたステータスを確認してください。

## ライセンス

MIT License
