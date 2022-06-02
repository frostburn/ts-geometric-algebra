# ts-geometric-algebra
Ts-geometric-algebra is a Clifford Algebra Generator for TypeScript and JavaScript. It generates Geometric Algebras of any signature and is mostly combatible with [Ganja.js](https://github.com/enkimute/ganja.js/). There is no operator overloading or inline functions yet due to lack of support in TypeScript.

(**Mathematically**, an algebra generated by ts-geometric-algebra is a graded exterior (Grassmann) algebra with a non-metric outer product, extended (Clifford) with geometric and contraction inner products, a Poincare duality operator and the main involutions and morphisms.)

(**Technically**, ts-geometric-algebra is a code generator producing classes that reificate algebraic literals.)

(**Practically**, ts-geometric-algebra enables algebraic operations over reals, complex numbers, dual numbers, hyperbolic numbers, vectors, spacetime events, quaternions, dual quaternions, biquaternions or any other Clifford Algebra with full type support from TypeScript.)

## Discourse and Discord

Visit [bivector.net](https://bivector.net) for our forum and chat - the perfect
place for questions and support.

### The Algebra Function

To create an Algebra, call the `Algebra` function specifying the metric
signature (number of positive, negative and zero dimensions). The result is
an ES6 class implementing the requested clifford algebra.

```typescript
function Algebra( p, q, r );
  // p    = number of positive dimensions.
  // q    = optional number of negative dimensions.
  // r    = optional number of zero dimensions.
```
An extended syntax is also available that allows you to further tweak the created Algebra.

```typescript
function Algebra( p, q, r, baseType, unroll )
  // p        = number of positive dimensions.
  // q        = optional number of negative dimensions.
  // r        = optional number of zero dimensions.
  // baseType = Float32Array or Float64Array
  // unroll   = (boolean) perform loop unrolling of performance critical methods
```
Here are some examples :

```typescript
// Basic
const Hyper   = Algebra(1);         // Hyperbolic numbers.
const Complex = Algebra(0, 1);      // Complex numbers.
const Dual    = Algebra(0, 0, 1);   // Dual numbers.
const H       = Algebra(0, 2);      // Quaternions.

// Clifford
const Cl2 = Algebra(2);             // Clifford algebra for 2D vector space.
const Cl3 = Algebra(3);             // Clifford algebra for 3D vector space.
const TimeSpace = Algebra(1, 3);    // Clifford algebra for timespace vectors.

// Geometric
const PGA2D = Algebra(2, 0, 1);     // Projective Euclidean 2D plane. (dual)
const PGA3D = Algebra(3, 0, 1);     // Projective Euclidean 3D space. (dual)
const CGA2D = Algebra(3, 1);        // Conformal 2D space.
const CGA3D = Algebra(4, 1);        // Conformal 3D space.

// High-Dimensional GA
const DCGA3D = Algebra(6, 2);       // Double Conformal 3D Space.
const TCGA3D = Algebra(9, 3);       // Triple Conformal 3D Space.
const DCGSTA = Algebra(4, 8);       // Double Conformal Geometric Space Time Algebra.
const QCGA   = Algebra(9, 6);       // Quadric Conformal Geometric Algebra.
```

You can now use these classes to generate algebraic elements. Those elements will have all of the
expected properties. (`norm`, blade access, `dot`, `wedge`, `mul`, `dual`, `inverse`, etc ...)

Unlike Ganja.js you must use them in a 'classic' programming style syntax like the example below.

```typescript
const Complex = Algebra(0, 1);       // Complex numbers.
const a = Complex.fromGanja([3, 2]); // 3 + 2i
const b = Complex.fromGanja([1, 4]); // 1 + 4i
return a.mul(b);                     // returns [-5, 14]
```
Altough not as pretty or fun as Ganja.js you have the full advantage of types and autocompletion.

### Methods

| Object oriented      | Function Oriented     | Ganja Equivalent | Explanation     |
|----------------------|-----------------------|------------------|-----------------|
| `x.equals(y)`        | `equals(x, y)`        | N/A              | Strict equality |
| `x.closeTo(y, tol?)` | `closeTo(x, y, tol?)` | N/A              | Equality within given tolerance |
| `x.hasNaN()`         | `hasNaN(x)`           | N/A              | Check for Not-a-Numbers |
| `x.hasInfinity()`    | `hasInfinity(x)`      | N/A              | Check for (negative) infinity |
| `x.isNil(tol?)`      | `isNil(x, tol?)`      | N/A              | Equal or close to zero |
| `x.isGrade(g, tol?)` | `isGrade(x, g, tol?)` | N/A              | Only has components of grade `g` |
| `x.s`                | N/A                   | `x.s`            | Scalar part (get/set) |
| `x.ps`               | N/A                   | N/A              | Pseudoscalar part (get/set) |
| `x.getAt(...idx)`    | N/A                   | N/A              | Metric-aware coefficient of the product of basis factors defined by `idx` |
| `x.setAt(...idx, a)` | N/A                   | N/A              | Set coefficient of product of basis factors `idx` as `a` |
| `x.norm()`           | `norm(x)`             | `x.Length`       | Conjugate norm (metric-aware) |
| `x.vnorm()`          | `vnorm(x)`            | `x.VLength`      | Vector norm (ignores metric) {`x.length` is array length} |
| `x.neg()`            | `neg(x)`              | N/A              | Negation (additive inverse) |
| `x.cwAbs()`          | `cwAbs(x)`            | N/A              | Component-wise absolute value |
| `x.involute()`       | `involute(x)`         | `x.Involute`     | Negation of basis factors |
| `x.rev()`            | `rev(x)`              | `x.Reverse`      | Reversal of basis factors {`x.reverse` is array reversal} |
| `x.conjugate()`      | `conjugate(x)`        | `x.Conjugate`    | Conjugation (combined involution and reversal) |
| `x.inverse()`        | `inverse(x)`          | `x.Inverse`      | Multiplicative inverse |
| `x.square()`         | `square(x)`           | `x.Mul(x)`       | Multiplicative squaring (optimized) |
| `x.normalize(a?)`    | `normalize(x, a?)`    | `x.Normalize`    | `x` with norm set to `a` (default 1) |
| `x.rotorNormalize()` | `rotorNormalize(x)`   | N/A              | Normalize rotor `x` |
| `x.sqrt()`           | `sqrt(x)`             | N/A              | Square root. Currently reliable only in dimensions < 2 |
| `x.rotorSqrt()`      | `rotorSqrt(x)`        | N/A              | Rotor square root. Available in certain metrics. |
| `x.exp()`            | `exp(x)`              | `x.Exp`          | Exponential function |
| `x.bivectorExp()`    | `bivectorExp(x)`      | `x.Exp`          | Bivector exponential function (optimized) |
| `x.log()`            | `log(x)`              | `x.Log`          | (Motor) Logarithm. Generic `exp` inverse available only in dimensions < 2 |
| `x.rotorLog()`       | `rotorLog(x)`         | `x.Log`          | Rotor Logarithm |
| `x.clone()`          | `clone(x)`            | `x.Scale(1)`     | Independent copy |
| `x.dual()`           | `dual(x)`             | `x.Dual` (*)     | Metric independent dual: `x.mul(x.dual()) = Cl.pseudoscalar()` |
| `x.undual()`         | `undual(x)`           | N/A              | Inverse of `x.dual()` |
| `x.scale(a)`         | `scale(x, a)`         | `x.Scale(a)`     | Scalar multiplication |
| `x.pow(n)`           | `pow(x, n)`           | `x.Pow(n)`       | Multiply `x` with itself `n` times |
| `x.applyWeights(ws)` | ...                   | ...              | ... |
| `x.negateGrades(...gs)` | ... | ... | ... |
| `x.add(y)` | ... | ... | ... |
| `x.sub(y)` | ... | ... | ... |
| `x.mul(y)` | ... | ... | ... |
| `x.rmul(y)` | ... | ... | ... |
| `x.div(y)` | ... | ... | ... |
| `x.ldiv(y)` | ... | ... | ... |
| `x.ldivs(y)` | ... | ... | ... |
| `x.wedge(y)` | ... | ... | ... |
| `x.rwedge(y)` | ... | ... | ... |
| `x.vee(y)` | ... | ... | ... |
| `x.rvee(y)` | ... | ... | ... |
| `x.rotorMean(y)` | ... | ... | ... |
| `x.contract(y, ctn)` | ... | ... | ... |
| `x.dot(y)` | ... | ... | ... |
| `x.dotL(y)` | ... | ... | ... |
| `x.dotR(y)` | ... | ... | ... |
| `x.star(y)` | ... | ... | ... |
| `x.imag()` | ... | ... | ... |
| `x.even()` | ... | ... | ... |
| `x.grade(1)` | ... | ... | ... |
| `x.grade(2)` | ... | ... | ... |
| `x.grade(n)` | ... | ... | ... |
| `x.vector()` | ... | ... | ... |
| `x.rotor()` | ... | ... | ... |
| `x.ganja()` | ... | ... | ... |
| `Cl.zero()` | ... | ... | ... |
| `Cl.scalar()` | ... | ... | ... |
| `Cl.pseudoscalar()` | ... | ... | ... |
| `Cl.basisVector(...idx)` | ... | ... | ... |
| `Cl.fromVector(vs, g)` | ... | ... | ... |
| `Cl.fromRotor(vs)` | ... | ... | ... |
| `Cl.fromGanja(vs)` | ... | ... | ... |
| `Cl.dimensions` | ... | ... | ... |
| `Cl.size` | ... | ... | ... |

(*) Only in degenerate metrics

### Dual Zoo
| Object oriented      | Function Oriented     | Ganja Equivalent    | Explanation     |
|----------------------|-----------------------|---------------------|-----------------|
| `x.podge()`          | `podge(x)`            | `x.Mul(ps)`         | Right-multiplication by `Cl.pseudoscalar()` |
| `x.unpodge()`        | `unpodge(x)`          | `x.Div(ps)`         | Right-division by `Cl.pseudoscalar()` |
| `x.podgeL()`         | `podgeL(x)`           | `ps.Mul(x)`         | Left-multiplication by `Cl.pseudoscalar()` |
| `x.unpodgeL()`       | `unpodgeL(x)`         | `ps.Inverse.Mul(x)` | Left-division by `Cl.pseudoscalar()` |
| `x.star()`           | `star(x)`             | N/A                 | Non-degenerate `x.podge()` |
| `x.unstar()`         | `unstar(x)`           | N/A                 | Inverse of `x.star()` |
| `x.starL(x)`         | `starL(x)`            | `x.Dual` (*)        | Non-degenerate `x.podgeL()` |
| `x.unstarL(x)`       | `unstarL(x)`          | N/A                 | Inverse of `x.starL()` |
| `x.hodge()`          | `hodge(x)`            | N/A                 | Hodge dual |
| `x.unhodge()`        | `unhodge(x)`          | N/A                 | Inverse of `x.hodge()` |
| `x.hodgeL()`         | `hodgeL(x)`           | N/A                 | Left Hodge dual |
| `x.unhodgeL()`       | `unhodgeL(x)`         | N/A                 | Inverse of `x.hodgeL()` |

(*) Only in non-degenerate metrics

#### Rotor operation availability
| Metric (pqr) | Available operations |
|--------------|----------------------|
| p+q+r <= 2   | (all)                |
| 400          | `rotorSqrt`, `rotorNormalize` |
| 310          | `rotorSqrt`, `rotorNormalize` |
| 301          | `rotorLog`, `rotorSqrt`, `rotorNormalize` |
| 410          | `rotorSqrt`, `rotorNormalize` |
