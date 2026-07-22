# minnesotanSpec

This page is retained for compatibility with older documentation links.

The current SVAR toolbox does not provide a `minnesotanSpec` class. Use
`minnesotabvarm` with `Type="normal"` or construct
`svar.minnesotanbvarm` directly.

## Current Alternatives

```matlab
PriorMdl = minnesotabvarm(size(Y,2),4,Y,Type="normal",lambda2=0.5);
PriorMdl = svar.minnesotanbvarm(size(Y,2),4,Y,lambda2=0.5);
```

## See Also

`minnesotabvarm`, `svar.minnesotanbvarm`, `svar.minnesotainwbvarm`
