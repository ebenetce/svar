function value = logDetPD(A)
%LOGDETPD Log-determinant of a symmetric positive-definite matrix.
R     = chol((A + A')/2);
value = 2*sum(log(diag(R)));
end