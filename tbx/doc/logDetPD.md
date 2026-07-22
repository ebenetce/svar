# svar.logDetPD

Log-determinant of a symmetric positive-definite matrix.

`svar.logDetPD` computes the log determinant of a positive-definite matrix
using a Cholesky factorization.

## Syntax

```matlab
value = svar.logDetPD(A)
```

## Description

`value = svar.logDetPD(A)` returns `log(det(A))` for a symmetric
positive-definite matrix `A`. The function symmetrizes `A` by averaging it with
its transpose before factorization.

## Input Arguments

`A` - Positive-definite matrix
: Numeric square matrix.

## Output Arguments

`value` - Log determinant
: Numeric scalar.

## Examples

### Compute a Log Determinant

```matlab
A = [2 0.5; 0.5 1];
value = svar.logDetPD(A);
```

## More About

### Numerical Method

The function uses `chol((A + A')/2)` and returns twice the sum of the logarithms
of the Cholesky diagonal entries. This avoids explicitly computing
`det(A)`.

## See Also

`chol`, `det`, `log`
