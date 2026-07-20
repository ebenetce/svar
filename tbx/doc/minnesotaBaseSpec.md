# minnesotaBaseSpec

Shared base class for Minnesota hyperparameter specifications.

`minnesotaBaseSpec` implements the common scalar-or-bounds behavior used by
the concrete Minnesota specification classes. This class is abstract and hidden;
construct specs with `minnesotaSpec` instead.

## Class Details

`minnesotaBaseSpec` is not intended for direct construction. It defines shared
properties, packing utilities, and validation used by `minnesotamniwSpec`,
`minnesotainwSpec`, and `minnesotanSpec`.

## Properties

`lambda1` - Overall tightness
: `0.2` (default) | positive scalar | two-element positive bounds |
  `hyperprior`.

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar | two-element positive bounds |
  `hyperprior`.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `[]` (default) | numeric row vector. An empty value is resolved to a vector
  of ones when a concrete model is built.

## Object Functions

`isFree`
: Return `true` when a named hyperparameter field is represented by bounds or a
  `hyperprior` object.

`freeFields`
: Return the names of all free hyperparameter fields.

`isResolved`
: Return `true` when all hyperparameter fields are scalar numeric values.

`pack`
: Return optimizer starting values, lower bounds, upper bounds, and field names
  for free fields.

`unpack`
: Write an optimizer vector back into selected fields and return the resolved
  spec.

`logHyperprior`
: Sum log-density terms for fields represented by `hyperprior` objects.

`hyperpriorFields`
: Return the names of fields that hold `hyperprior` objects.

## More About

### Scalar-or-Bounds Convention

Minnesota spec fields use the value itself to determine whether a field is
fixed or free:

* A scalar value fixes the field.
* A two-element vector `[lower upper]` makes the field free with flat bounds.
* A `hyperprior` object makes the field free and adds a log-prior term.

This convention keeps the specification object as the single source of truth
for optimization.

### Build Contract

Concrete `build` methods call `assertResolvedForBuild`. A spec can be built only
after every free field has been replaced by a scalar candidate value.

## See Also

`minnesotaSpec`, `minnesotamniwSpec`, `minnesotainwSpec`, `minnesotanSpec`,
`hyperprior`

