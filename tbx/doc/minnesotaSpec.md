# minnesotaSpec

This page is retained for compatibility with older documentation links.

The current SVAR toolbox does not provide a `minnesotaSpec` function. Use
`minnesotabvarm` to create Minnesota priors directly, or use `glp` to tune the
conjugate Minnesota family by marginal likelihood.

## Current Alternatives

```matlab
PriorMdl = minnesotabvarm(size(Y,2),4,Y);
TunedMdl = glp(size(Y,2),4,Y,Psi="conditional");
```

Use `Type="inw"` or `Type="normal"` with `minnesotabvarm` to create the
semiconjugate or fixed-Sigma families.

## See Also

`minnesotabvarm`, `glp`, `svar.minnesotamniwbvarm`,
`svar.minnesotainwbvarm`, `svar.minnesotanbvarm`
