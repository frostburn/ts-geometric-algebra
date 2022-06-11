import {eigenValues, sinc, sinch} from './utils';
import {type ElementBaseType, type AlgebraElement} from './element';
import {pqrMixin} from './pqr';
import {linSolve} from './element';

export * from './element';

export function vLinSolve(x: number[], basis: number[][], threshold = 1e-6) {
  const Grassmann = Algebra(0, 0, x.length);
  return linSolve(
    Grassmann.fromVector(x),
    basis.map(b => Grassmann.fromVector(b)),
    threshold
  );
}

// https://stackoverflow.com/a/43122214
function bitCount(n: number) {
  n = n - ((n >> 1) & 0x55555555);
  n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
  return (((n + (n >> 4)) & 0xf0f0f0f) * 0x1010101) >> 24;
}

// Split number n into powers of two such that n = f1 ^ f2 ^ ... ^ fk
function bitSplit(n: number) {
  if (n < 0) {
    throw new Error('Cannot split negative bitsets');
  }
  const result = [];
  let p = 1;
  while (n) {
    if (n & 1) {
      result.push(p);
    }
    n >>= 1;
    p <<= 1;
  }
  return result;
}

// https://en.wikipedia.org/wiki/Gnome_sort
function sortSign(sequence: number[]) {
  let sign = 1;
  let pos = 0;
  while (pos < sequence.length) {
    if (pos === 0 || sequence[pos] >= sequence[pos - 1]) {
      pos++;
    } else {
      const temp = sequence[pos];
      sequence[pos] = sequence[pos - 1];
      sequence[pos - 1] = temp;
      sign = -sign;
      pos--;
    }
  }
  return sign;
}

// Contraction criteria
function symmetric(r: number, s: number) {
  return Math.abs(r - s);
}

function left(r: number, s: number) {
  return s - r;
}

function right(r: number, s: number) {
  return r - s;
}

// eslint-disable-next-line @typescript-eslint/no-unused-vars
function nil(r: number, s: number) {
  return 0;
}

function reduceIndices(indices: number[]) {
  return indices.map(i => 1 << i).reduce((a, b) => a ^ b, 0);
}

const MAX_DIMENSIONS = 36;

export function Algebra(
  p: number,
  q = 0,
  r = 0,
  baseType: typeof ElementBaseType = Float32Array,
  unroll = true
): typeof AlgebraElement {
  const metric: number[] = [];
  for (let i = 0; i < r; ++i) {
    metric.push(0);
  }
  for (let i = 0; i < p; ++i) {
    metric.push(1);
  }
  for (let i = 0; i < q; ++i) {
    metric.push(-1);
  }
  const dimensions = p + q + r;
  const size = 1 << dimensions;
  const indexMask = size - 1;

  // Helper class used in factorization for easier inverses
  let Euclidized: typeof AlgebraElement;

  function getEuclidized() {
    if (Euclidized === undefined) {
      Euclidized = Algebra(dimensions, 0, 0, baseType, unroll);
    }
    return Euclidized;
  }

  if (dimensions > MAX_DIMENSIONS) {
    throw new Error(`Maximum total number of dimensions is ${MAX_DIMENSIONS}`);
  }

  // Geometric product between basis vectors
  function basisIndexMul(...indices: number[]) {
    // sign incorporates Ex * Ey metric
    const sign = sortSign(indices);
    // weight incorporates Ex * Ex metric
    let weight = 1;
    let i = 1;
    while (i < indices.length) {
      if (indices[i] === indices[i - 1]) {
        weight *= metric[indices.splice(i, 1)[0]];
      } else {
        i++;
      }
    }
    return sign * weight;
  }

  // Geometric product between bundles of basis vectors using bit field indices
  function basisMul(...bitFieldIndices: number[]) {
    const indices: number[] = [];

    bitFieldIndices.forEach(bitFieldIndex => {
      for (let i = 0; i < dimensions; ++i) {
        const p = 1 << i;
        if (p & bitFieldIndex) {
          indices.push(i);
        }
      }
    });
    return basisIndexMul(...indices);
  }

  // This could be turned into a bit array if memory becomes an issue
  const mulTable: number[][] = [];
  for (let i = 0; i < size; ++i) {
    const row: number[] = [];
    for (let j = 0; j < size; ++j) {
      row.push(basisMul(i, j));
    }
    mulTable.push(row);
  }

  // Mapping from bit-field indices to ganja.js lexicographic order
  const indexString: [number, string][] = [];
  for (let i = 0; i < size; ++i) {
    let str = '';
    for (let j = 0; j < dimensions; ++j) {
      if (i & (1 << j)) {
        str += j.toString(MAX_DIMENSIONS);
      }
    }
    indexString.push([i, str]);
  }
  function cmp(a: [any, string], b: [any, string]) {
    if (a[1].length < b[1].length) {
      return -1;
    }
    if (a[1].length > b[1].length) {
      return 1;
    }
    if (a[1] < b[1]) {
      return -1;
    }
    if (a[1] > b[1]) {
      return 1;
    }
    return 0;
  }
  indexString.sort(cmp);

  class AlgebraClass extends baseType {
    constructor(values?: Iterable<number>) {
      if (values === undefined) {
        super(size);
      } else {
        super(values);
      }
    }

    // This is a hack to get around TypeScript's lack of abstract static methods
    get cls() {
      return AlgebraClass;
    }

    empty() {
      return new this.cls();
    }

    zeroed() {
      return this.cls.zero();
    }

    equals(other: AlgebraElement) {
      for (let i = 0; i < this.length; ++i) {
        if (this[i] !== other[i]) {
          return false;
        }
      }
      return true;
    }

    closeTo(other: AlgebraElement, tolerance = 1e-4) {
      for (let i = 0; i < this.length; ++i) {
        const error = Math.abs(this[i] - other[i]);
        if (error > tolerance || isNaN(error)) {
          return false;
        }
      }
      return true;
    }

    hasNaN() {
      for (let i = 0; i < this.length; ++i) {
        if (isNaN(this[i])) {
          return true;
        }
      }
      return false;
    }

    hasInfinity() {
      for (let i = 0; i < this.length; ++i) {
        if (Math.abs(this[i]) === Infinity) {
          return true;
        }
      }
      return false;
    }

    isNil(tolerance = 0) {
      for (let i = 0; i < this.length; ++i) {
        if (Math.abs(this[i]) > tolerance) {
          return false;
        }
      }
      return true;
    }

    isGrade(grade: number, tolerance = 0) {
      for (let i = 0; i < this.length; ++i) {
        if (bitCount(i) !== grade && Math.abs(this[i]) > tolerance) {
          return false;
        }
      }
      return true;
    }

    norm() {
      return Math.sqrt(Math.abs(this.mul(this.conjugate()).s));
    }

    vnorm() {
      let result = 0;
      for (let i = 0; i < this.length; ++i) {
        result += this[i] * this[i];
      }
      return Math.sqrt(result);
    }

    neg(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = -this[i];
      }
      return result;
    }

    cwAbs(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = Math.abs(this[i]);
      }
      return result;
    }

    rev(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = bitCount(i) & 2 ? -this[i] : this[i];
      }
      return result;
    }

    involute(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = bitCount(i) & 1 ? -this[i] : this[i];
      }
      return result;
    }

    conjugate(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = (bitCount(i) + 1) & 2 ? -this[i] : this[i];
      }
      return result;
    }

    // === Dual Zoo ===

    // For all Ex = AlgebraClass.basisBlade(...x)
    // Ex.mul(Ex.dual()) === AlgebraClass.pseudoscalar()
    dual(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] = this[i] * mulTable[i][dualIndex];
      }
      return result;
    }
    // a.vee(b) === undual(b.dual().wedge(a.dual()))
    undual(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] = this[i] * mulTable[dualIndex][i];
      }
      return result;
    }

    podge(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[indexMask ^ i] = this[i] * mulTable[i][indexMask];
      }
      return result;
    }
    unpodge(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[indexMask ^ i] = this[i] / mulTable[indexMask ^ i][indexMask];
      }
      return result;
    }
    podgeL(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[indexMask ^ i] = this[i] * mulTable[indexMask][i];
      }
      return result;
    }
    unpodgeL(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[indexMask ^ i] = this[i] / mulTable[indexMask][indexMask ^ i];
      }
      return result;
    }

    // See star() overload for forward implementation
    unstar(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[indexMask ^ i] =
          this[i] * (mulTable[indexMask ^ i][indexMask] || 1);
      }
      return result;
    }
    starL(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[indexMask ^ i] = this[i] * (mulTable[indexMask][i] || 1);
      }
      return result;
    }
    unstarL(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[indexMask ^ i] =
          this[i] * (mulTable[indexMask][indexMask ^ i] || 1);
      }
      return result;
    }

    hodge(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] =
          this[i] * (mulTable[i][indexMask] || mulTable[i][dualIndex]);
      }
      return result;
    }
    unhodge(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] =
          this[i] * (mulTable[dualIndex][indexMask] || mulTable[dualIndex][i]);
      }
      return result;
    }
    hodgeL(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] =
          this[i] * (mulTable[indexMask][i] || mulTable[dualIndex][i]);
      }
      return result;
    }
    unhodgeL(): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] =
          this[i] * (mulTable[indexMask][dualIndex] || mulTable[i][dualIndex]);
      }
      return result;
    }

    normalize(newNorm = 1): AlgebraElement {
      return this.scale(newNorm / this.norm());
    }

    rotorNormalize(): AlgebraElement {
      if (dimensions <= 2) {
        return this.normalize();
      }
      throw new Error('Do not know how to normalize a rotor in this algebra');
    }

    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    sqrt(forceBabylon = false, numIter = 16): AlgebraElement {
      // Not quaranteed to converge. Barely better than nothing.
      // https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Babylonian_method
      const root = this.plus(1);
      for (let i = 1; i < numIter; ++i) {
        root.accumulate(this.div(root)).rescale(0.5);
      }
      return root;
    }

    rotorSqrt(): AlgebraElement {
      const root = this.plus(1);
      return root.rotorNormalize();
    }

    exp(forceTaylor = false, numTaylorTerms = 32): AlgebraElement {
      if (!forceTaylor) {
        // Closed form exp
        const grade2 = this.imag();
        if (grade2.isGrade(2)) {
          return grade2.split().reduce((total, simple) => {
            const square = simple.square().s;
            const len = Math.sqrt(Math.abs(square));
            if (len <= 1e-12) {
              simple.s += 1;
            } else if (square < 0) {
              simple = simple.scale(sinc(len));
              simple.s += Math.cos(len);
            } else {
              simple = simple.scale(sinch(len));
              simple.s += Math.cosh(len);
            }
            return total.mul(simple);
          }, this.cls.scalar(Math.exp(this.s)));
        }
      }

      // No specific implementation found, but we can still
      // use the fact that the scalar commutes with everything
      let maybeImag: AlgebraElement;
      if (forceTaylor) {
        maybeImag = this;
      } else {
        maybeImag = this.imag();
        // In 3D the pseudoscalar always commutes
        if (dimensions === 3) {
          maybeImag.ps = 0;
        }
      }

      // Taylor series
      const result = this.cls.scalar();
      let term = this.cls.scalar();
      for (let i = 1; i < numTaylorTerms; ++i) {
        term = term.mul(maybeImag.scale(1 / i));
        result.accumulate(term);
      }

      if (forceTaylor) {
        return result;
      }
      if (dimensions === 3) {
        const factor = this.cls.zero();
        factor.s = factor.ps = Math.exp(this.s);
        if (mulTable[indexMask][indexMask] > 0) {
          factor.s *= Math.cosh(this.ps);
          factor.ps *= Math.sinh(this.ps);
        } else if (mulTable[indexMask][indexMask] < 0) {
          factor.s *= Math.cos(this.ps);
          factor.ps *= Math.sin(this.ps);
        } else {
          factor.ps *= this.ps;
        }
        return result.mul(factor);
      }
      return result.rescale(Math.exp(this.s));
    }

    bivectorExp() {
      return this.exp();
    }

    log(forceProduct = false, numProductTerms = 20): AlgebraElement {
      if (!forceProduct && this.isGrade(2)) {
        const sum = this.zeroed();
        this.motorFactorize().forEach(bi => {
          const [ci, si] = [bi.s, bi.grade(2)];
          const square = si.square().s;
          const len = Math.sqrt(Math.abs(square));
          if (Math.abs(square) < 1e-5) sum.accumulate(si);
          if (square < 0) sum.accumulate(si.scale(Math.acos(ci) / len));
          sum.accumulate(si.scale(Math.acosh(ci) / len));
        });
        return sum;
      }
      // Simply assumed to work for general multi vectors, someone should prove this
      // https://www.emis.de/journals/HOA/IJMMS/2004/65-683653.pdf
      let result = this.plus(-1);
      let term = this.clone();
      for (let i = 0; i < numProductTerms; ++i) {
        term = term.sqrt();
        result = result.mul(term.plus(1).inverse().rescale(2));
      }
      return result;
    }

    rotorLog() {
      return this.log();
    }

    clone(): AlgebraElement {
      return new this.cls(this);
    }

    negateGrades(...grades: number[]): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = grades.includes(bitCount(i)) ? -this[i] : this[i];
      }
      return result;
    }

    adjugate(): AlgebraElement {
      switch (dimensions) {
        case 0:
          return this.clone();
        case 1:
          return this.involute();
        case 2:
          return this.conjugate();
        case 3:
          return this.rev().mul(this.involute()).mul(this.conjugate());
        case 4:
          const conjugate = this.conjugate();
          return conjugate.mul(this.mul(conjugate).negateGrades(3, 4));
        case 5:
          const civ = this.conjugate().mul(this.involute()).mul(this.rev());
          return civ.mul(this.mul(civ).negateGrades(1, 4));
        default:
          // Shirokov-style adjugate
          const N = 1 << (((dimensions + 1) / 2) | 0);
          let Uk = this.clone();
          let adjU = this.clone();
          for (let k = 1; k < N; ++k) {
            let s = Uk.s * N;
            if (s % k !== 0) {
              s *= k;
              Uk.rescale(k);
              adjU.rescale(k);
            }
            adjU = Uk.plus(-s / k);
            Uk = this.mul(adjU);
          }
          return adjU;
      }
    }

    inverse(): AlgebraElement {
      // Matrix-free inverses up to 5D
      // http://repository.essex.ac.uk/17282/1/TechReport_CES-534.pdf
      switch (dimensions) {
        // cases 0, 1 and 2 completely covered by pqr.ts
        case 3:
          const reverse = this.rev();
          const involute = this.involute();
          const conjugate = this.conjugate();
          return reverse
            .mul(involute)
            .mul(conjugate)
            .scale(1 / this.mul(conjugate).mul(involute).mul(reverse).s);
        case 4:
          const conjugate4 = this.conjugate();
          const modulus = this.mul(conjugate4);
          const n34 = modulus.negateGrades(3, 4);
          return conjugate4.mul(n34).scale(1 / modulus.mul(n34).s);
        case 5:
          const civ = this.conjugate().mul(this.involute()).mul(this.rev());
          const tciv = this.mul(civ);
          const tciv14 = tciv.negateGrades(1, 4);
          return civ.mul(tciv14).scale(1 / tciv.mul(tciv14).s);
        default:
          // Shirokov inverse
          const N = 1 << (((dimensions + 1) / 2) | 0);
          let Uk = this.clone();
          let adjU: AlgebraElement;
          for (let k = 1; k < N; ++k) {
            adjU = Uk.plus(-(N / k) * Uk.s);
            Uk = this.mul(adjU);
          }
          return Uk.s === 0 ? this.zeroed() : adjU!.scale(1 / Uk.s);
      }
    }

    square(): AlgebraElement {
      return this.mul(this);
    }

    scale(scalar: number): AlgebraElement {
      return this.clone().rescale(scalar);
    }

    rescale(scalar: number): AlgebraElement {
      for (let i = 0; i < this.length; ++i) {
        this[i] *= scalar;
      }
      return this;
    }

    pow(power: number, splitStages = 8): AlgebraElement {
      if (power !== Math.round(power)) {
        if (dimensions === 0) {
          return this.cls.scalar(Math.pow(this.s, power));
        } else if (dimensions === 1) {
          return this.log().scale(power).exp();
        } else if (p === 0 && q === 2 && r === 0) {
          return this.log().scale(power).exp();
        }
        let epsilon = this.sqrt();
        power *= 2;
        let stages = splitStages;
        while (stages > 0 && power !== Math.round(power)) {
          epsilon = epsilon.sqrt();
          power *= 2;
          stages--;
        }
        return epsilon.pow(power);
      }
      if (power === 0) {
        return this.cls.scalar();
      }
      if (power === 1) {
        return this.clone();
      }
      if (power === 2) {
        return this.square();
      }
      if (power > 0) {
        let result = this.cls.scalar();
        let powerOfTwo = this.clone();
        while (power) {
          if (power & 1) {
            result = result.mul(powerOfTwo);
          }
          powerOfTwo = powerOfTwo.square();
          power >>= 1;
        }
        return result;
      }
      return this.inverse().pow(-power);
    }

    applyWeights(weights: number[]) {
      const result = this.clone();
      for (let i = 0; i < this.length; ++i) {
        for (let j = 0; j < weights.length; ++j) {
          if (i & (1 << j)) {
            result[i] *= weights[j];
          }
        }
      }
      return result;
    }

    add(other: AlgebraElement): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i] + other[i];
      }
      return result;
    }

    sub(other: AlgebraElement): AlgebraElement {
      const result = this.empty();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i] - other[i];
      }
      return result;
    }

    mul(other: AlgebraElement): AlgebraElement {
      const result = this.zeroed();
      for (let i = 0; i < this.length; ++i) {
        if (!this[i]) {
          continue;
        }
        for (let j = 0; j < other.length; ++j) {
          result[i ^ j] += this[i] * other[j] * mulTable[i][j];
        }
      }
      return result;
    }

    lmul(other: AlgebraElement): AlgebraElement {
      const result = this.zeroed();
      for (let i = 0; i < this.length; ++i) {
        if (!this[i]) {
          continue;
        }
        for (let j = 0; j < other.length; ++j) {
          result[i ^ j] += this[i] * other[j] * mulTable[j][i];
        }
      }
      return result;
    }

    div(other: AlgebraElement): AlgebraElement {
      return this.mul(other.inverse());
    }

    ldiv(other: AlgebraElement): AlgebraElement {
      return other.inverse().mul(this);
    }

    ldivs(other: AlgebraElement): AlgebraElement {
      return this.inverse().mul(other);
    }

    wedge(other: AlgebraElement): AlgebraElement {
      const result = this.zeroed();
      for (let i = 0; i < this.length; ++i) {
        if (!this[i]) {
          continue;
        }
        for (let j = 0; j < other.length; ++j) {
          if (!(i & j)) {
            result[i ^ j] += this[i] * other[j] * mulTable[i][j];
          }
        }
      }
      return result;
    }

    lwedge(other: AlgebraElement): AlgebraElement {
      const result = this.zeroed();
      for (let i = 0; i < this.length; ++i) {
        if (!this[i]) {
          continue;
        }
        for (let j = 0; j < other.length; ++j) {
          if (!(i & j)) {
            result[i ^ j] += this[i] * other[j] * mulTable[j][i];
          }
        }
      }
      return result;
    }

    vee(other: AlgebraElement): AlgebraElement {
      return other.dual().wedge(this.dual()).undual();
    }

    lvee(other: AlgebraElement): AlgebraElement {
      return other.dual().lwedge(this.dual()).undual();
    }

    delta(other: AlgebraElement, threshold = 0): AlgebraElement {
      const product = this.mul(other);
      const maxGrade = Math.max(0, ...product.grades(threshold));
      return product.grade(maxGrade);
    }

    rotorMean(other: AlgebraElement) {
      return this.add(other).rotorNormalize();
    }

    contract(
      other: AlgebraElement,
      criterion: (r: number, s: number) => number
    ): AlgebraElement {
      const result = this.zeroed();
      for (let i = 0; i < this.length; ++i) {
        if (!this[i]) {
          continue;
        }
        const gradeI = bitCount(i);
        for (let j = 0; j < other.length; ++j) {
          const gradeJ = bitCount(j);
          const target = i ^ j;
          const gradeTarget = bitCount(target);
          if (gradeTarget === criterion(gradeI, gradeJ)) {
            result[target] += this[i] * other[j] * mulTable[i][j];
          }
        }
      }
      return result;
    }

    dot(other: AlgebraElement) {
      return this.contract(other, symmetric);
    }

    dotL(other: AlgebraElement) {
      return this.contract(other, left);
    }

    ldotL(other: AlgebraElement) {
      return other.contract(this, left);
    }

    dotR(other: AlgebraElement) {
      return this.contract(other, right);
    }

    ldotR(other: AlgebraElement) {
      return other.contract(this, right);
    }

    dotS(other: AlgebraElement) {
      return this.contract(other, nil);
    }

    // Scalar part
    get s(): number {
      return this[0];
    }
    set s(value: number) {
      this[0] = value;
    }

    // Pseudoscalar part
    get ps(): number {
      return this[this.length - 1];
    }
    set ps(value: number) {
      this[this.length - 1] = value;
    }

    getAt(...indices: number[]): number {
      return this[reduceIndices(indices)] * basisIndexMul(...indices);
    }

    setAt(...indicesAndValue: number[]): this {
      const indices = indicesAndValue.slice(0, -1);
      const value = indicesAndValue[indicesAndValue.length - 1];
      this[reduceIndices(indices)] = value / basisIndexMul(...indices);
      return this;
    }

    imag() {
      const result = this.clone();
      result.s = 0;
      return result;
    }

    even() {
      const result = this.zeroed();
      for (let i = 0; i < this.length; ++i) {
        if (!(bitCount(i) & 1)) {
          result[i] = this[i];
        }
      }
      return result;
    }

    grade(grade: number) {
      const result = this.zeroed();
      for (let i = 0; i < this.length; ++i) {
        if (bitCount(i) === grade) {
          result[i] = this[i];
        }
      }
      return result;
    }

    vector(grade = 1) {
      const result = [];
      for (let i = 0; i < this.length; ++i) {
        if (indexString[i][1].length === grade) {
          result.push(this[indexString[i][0]]);
        }
      }
      return new baseType(result);
    }

    rotor() {
      const result = [];
      for (let i = 0; i < this.length; ++i) {
        if (indexString[i][1].length % 2 === 0) {
          result.push(this[indexString[i][0]]);
        }
      }
      return new baseType(result);
    }

    // Ganja.js compatible representation
    ganja() {
      return new baseType(indexString.map(g => this[g[0]]));
    }

    invScale(other: AlgebraElement, threshold = 0) {
      let result = NaN;
      for (let i = 0; i < this.length; ++i) {
        if (Math.abs(other[i]) > threshold) {
          if (isNaN(result)) {
            result = this[i] / other[i];
          } else {
            if (Math.abs(1 - this[i] / other[i] / result) > threshold) {
              return NaN;
            }
          }
        }
      }
      return result;
    }

    grades(threshold = 0) {
      const grades = new Set<number>();
      for (let i = 0; i < this.length; ++i) {
        if (Math.abs(this[i]) > threshold) {
          grades.add(bitCount(i));
        }
      }
      const result = [...grades.values()];
      result.sort((a, b) => a - b);
      return result;
    }

    plus(scalar: number) {
      const result = this.clone();
      result.s += scalar;
      return result;
    }

    accumulate(other: AlgebraElement) {
      for (let i = 0; i < this.length; ++i) {
        this[i] += other[i];
      }
      return this;
    }

    bladeFactorize(): [AlgebraElement[], number] {
      const Euclid = getEuclidized();
      let euclid = new Euclid(this);
      const norm = euclid.norm();
      euclid.rescale(1 / norm);

      let maxIndex = 0;
      let maxCoordinate = euclid[0];
      for (let i = 1; i < euclid.length; ++i) {
        if (euclid[i] > maxCoordinate) {
          maxIndex = i;
          maxCoordinate = euclid[i];
        }
      }
      const indices = bitSplit(maxIndex);
      const factors = [];
      for (let i = 0; i < indices.length - 1; ++i) {
        const basisBlade = Euclid.zero();
        basisBlade[indices[i]] = 1;
        const factor = basisBlade.dotL(euclid.inverse()).dotL(euclid);
        factor.rescale(1 / factor.norm());
        factors.push(new this.cls(factor));
        euclid = factor.inverse().dotL(euclid);
      }
      euclid.rescale(1 / euclid.norm());
      factors.push(new this.cls(euclid));

      return [factors, norm];
    }

    // Bivector split - we handle all real cases, still have to add the complex cases for those exception scenarios.
    split(iter = 50) {
      const TWT = this.wedge(this);
      if (TWT.vnorm() < 1e-5) return [this.clone()]; // bivector was simple.
      const k = Math.floor(dimensions / 2);
      let B = this.clone();
      let m = 1;
      let Wi: AlgebraElement[] = [];
      for (let i = 0; i < k; ++i) {
        m = m * (i + 1);
        Wi.push(B.scale(1 / m));
        B = B.wedge(this);
      }
      let eigen;
      if (k < 3) {
        // The quadratic case is easy to solve. (for spaces <6D)
        const TDT = this.dot(this).s;
        const D = 0.5 * Math.sqrt(TDT * TDT - TWT.square().s);
        eigen = [0.5 * TDT + D, 0.5 * TDT - D].sort(
          (a, b) => Math.abs(a) - Math.abs(b)
        );
      } else {
        // For >6D, closed form solutions of the characteristic polyn. are impossible, use eigenvalues of companion matrix.
        const Wis = Wi.map((W, i) => W.square().s * (-1) ** (k - i + (k % 2)));
        const matrix: number[][] = [];
        for (let i = 0; i < k; ++i) {
          const row: number[] = [];
          for (let j = 0; j < k; ++j) {
            if (j === k - 1) {
              row.push(Wis[k - i - 1]);
            } else if (i - 1 === j) {
              row.push(1);
            } else {
              row.push(0);
            }
          }
          matrix.push(row);
        }
        eigen = eigenValues(matrix, iter).sort(
          (a, b) => Math.abs(a) - Math.abs(b)
        );
      }
      Wi = [this.cls.scalar(), ...Wi, this.zeroed()];
      const sum = this.zeroed();
      const k2 = Math.floor(k / 2);
      const res: AlgebraElement[] = eigen.slice(1).map(v => {
        const N = this.zeroed();
        const DN = this.zeroed();
        for (let i = 0; i <= k2; ++i) {
          N.accumulate(Wi[2 * i + 1].scale(v ** (k2 - i)));
          DN.accumulate(Wi[2 * i].scale(v ** (k2 - i)));
        }
        if (DN.vnorm() === 0) return this.zeroed();
        const ret = N.div(DN);
        sum.accumulate(ret);
        return ret;
      });
      return [this.sub(sum), ...res]; // Smallest eigvalue becomes B-rest
    }

    // Factorize a motor
    motorFactorize(iter = 50) {
      const S = this.grade(2).split(iter);
      const R = S.slice(0, S.length - 1).map(Mi => {
        Mi.s += this.s;
        const scale = Math.sqrt(Mi.rev().mul(Mi).s);
        return Mi.scale(1 / scale);
      });
      R.push(
        R.reduce((tot, fact) => tot.mul(fact.rev()), this.cls.scalar()).mul(
          this
        )
      );
      return R;
    }

    meetJoin(
      other: AlgebraElement,
      threshold = 0
    ): [AlgebraElement, AlgebraElement] {
      const grades = this.grades(threshold);
      const otherGrades = other.grades(threshold);
      if (grades.length !== 1 || otherGrades.length !== 1) {
        throw new Error('Inputs must be blades');
      }
      const Euclid = getEuclidized();

      let a: AlgebraElement;
      let b: AlgebraElement;
      let aGrade: number;
      let bGrade: number;
      if (grades[0] <= otherGrades[0]) {
        a = new Euclid(this);
        b = new Euclid(other);
        aGrade = grades[0];
        bGrade = otherGrades[0];
      } else {
        a = new Euclid(other);
        b = new Euclid(this);
        aGrade = otherGrades[0];
        bGrade = grades[0];
      }
      const d = a.delta(b, threshold);
      const dGrade = d.grades(threshold)[0];
      const joinGrade = (aGrade + bGrade + dGrade) / 2;
      const meetGrade = (aGrade + bGrade - dGrade) / 2;

      const factors = d.dual().bladeFactorize()[0];

      let meet = Euclid.scalar();
      let join = Euclid.pseudoscalar();
      let mg = 0;
      let jg = dimensions;

      const aInverse = a.inverse();
      for (let i = 0; i < factors.length; ++i) {
        const factor = factors[i];
        const projection = factor.dotL(aInverse).dotL(a);
        const rejection = factor.sub(projection);
        if (!projection.isNil(threshold)) {
          meet = meet.wedge(projection);
          mg++;
          if (mg === meetGrade) {
            join = a.wedge(meet.inverse().dotL(b));
            break;
          }
        }
        if (!rejection.isNil(threshold)) {
          join = rejection.dotL(join);
          jg--;
          if (jg === joinGrade) {
            meet = b.dotL(join.inverse()).dotL(a);
            break;
          }
        }
      }

      return [
        new this.cls(meet).grade(meetGrade),
        new this.cls(join).grade(joinGrade),
      ];
    }

    star(): AlgebraElement; // Dischord dual
    star(other: AlgebraElement): number; // Scalar product
    star(maybeOther?: AlgebraElement) {
      if (maybeOther === undefined) {
        const result = this.empty();
        for (let i = 0; i < this.length; ++i) {
          result[indexMask ^ i] = this[i] * (mulTable[i][indexMask] || 1);
        }
        return result as AlgebraElement;
      }
      let result = 0;
      for (let i = 0; i < this.length; ++i) {
        result += this[i] * maybeOther[i] * mulTable[i][i];
      }
      return result;
    }

    static zero(): AlgebraElement {
      return new AlgebraClass().fill(0);
    }

    static scalar(magnitude = 1): AlgebraElement {
      const result = AlgebraClass.zero();
      result[0] = magnitude;
      return result;
    }

    static pseudoscalar(magnitude = 1): AlgebraElement {
      const result = AlgebraClass.zero();
      result[size - 1] = magnitude;
      return result;
    }

    static basisBlade(...indices: number[]): AlgebraElement {
      const result = AlgebraClass.zero();
      result[reduceIndices(indices)] = 1 / basisIndexMul(...indices);
      return result;
    }

    static fromVector(values: Iterable<number>, grade = 1) {
      const result = AlgebraClass.zero();
      let i = 0;
      for (const component of values) {
        while (i < size) {
          if (indexString[i][1].length === grade) {
            result[indexString[i][0]] = component;
            i++;
            break;
          }
          i++;
        }
      }
      return result;
    }

    static fromRotor(values: Iterable<number>) {
      const result = AlgebraClass.zero();
      let i = 0;
      for (const component of values) {
        while (i < size) {
          if (indexString[i][1].length % 2 === 0) {
            result[indexString[i][0]] = component;
            i++;
            break;
          }
          i++;
        }
      }
      return result;
    }

    static fromGanja(values: Iterable<number>) {
      const result = new AlgebraClass();
      let index = 0;
      for (const component of values) {
        result[indexString[index++][0]] = component;
      }
      return result;
    }

    static get size() {
      return size;
    }

    static get dimensions() {
      return dimensions;
    }
  }

  const Result = pqrMixin(p, q, r, AlgebraClass);

  if (p === dimensions) {
    Euclidized = Result;
  }

  if (!unroll) {
    return Result;
  }

  // === Replace generic code with optimized unrolled versions ===

  let addInner = '';
  let subInner = '';
  for (let i = 0; i < size; ++i) {
    addInner += `res[${i}]=t[${i}]+o[${i}];`;
    subInner += `res[${i}]=t[${i}]-o[${i}];`;
  }

  const mulLines: string[] = [];
  for (let i = 0; i < size; ++i) {
    mulLines.push(`res[${i}]=`);
  }
  const wedgeLines = [...mulLines];
  const veeLines = [...mulLines];
  const dotLines = [...mulLines];
  const dotLeftLines = [...mulLines];

  for (let i = 0; i < size; ++i) {
    for (let j = 0; j < size; ++j) {
      if (mulTable[i][j] > 0) {
        mulLines[i ^ j] += `+t[${i}]*o[${j}]`;
        if (!(i & j)) {
          wedgeLines[i ^ j] += `+t[${i}]*o[${j}]`;
          veeLines[indexMask ^ i ^ j] += `+t[${indexMask ^ i}]*o[${
            indexMask ^ j
          }]`;
        }
      } else if (mulTable[i][j] < 0) {
        mulLines[i ^ j] += `-t[${i}]*o[${j}]`;
        if (!(i & j)) {
          wedgeLines[i ^ j] += `-t[${i}]*o[${j}]`;
          veeLines[indexMask ^ i ^ j] += `-t[${indexMask ^ i}]*o[${
            indexMask ^ j
          }]`;
        }
      }
    }
  }

  for (let i = 0; i < size; ++i) {
    const gradeI = bitCount(i);
    for (let j = 0; j < size; ++j) {
      const gradeJ = bitCount(j);
      const target = i ^ j;
      const gradeTarget = bitCount(target);
      if (mulTable[i][j] > 0) {
        if (gradeTarget === symmetric(gradeI, gradeJ)) {
          dotLines[target] += `+t[${i}]*o[${j}]`;
        }
        if (gradeTarget === left(gradeI, gradeJ)) {
          dotLeftLines[target] += `+t[${i}]*o[${j}]`;
        }
      } else if (mulTable[i][j] < 0) {
        if (gradeTarget === symmetric(gradeI, gradeJ)) {
          dotLines[target] += `-t[${i}]*o[${j}]`;
        }
        if (gradeTarget === left(gradeI, gradeJ)) {
          dotLeftLines[target] += `-t[${i}]*o[${j}]`;
        }
      }
    }
  }

  const squareTerms: Map<string, number>[] = [];
  for (let i = 0; i < size; ++i) {
    squareTerms.push(new Map());
  }

  for (let i = 0; i < size; ++i) {
    for (let j = 0; j < size; ++j) {
      const terms = squareTerms[i ^ j];
      const key = i < j ? `t[${i}]*t[${j}]` : `t[${j}]*t[${i}]`;
      if (terms.has(key)) {
        terms.set(key, terms.get(key)! + mulTable[i][j]);
      } else {
        terms.set(key, mulTable[i][j]);
      }
    }
  }
  const squareLines: string[] = [];
  for (let i = 0; i < size; ++i) {
    let ones = '';
    let twos = '';
    for (const [key, value] of squareTerms[i].entries()) {
      switch (value) {
        case 0:
          break;
        case 1:
          ones += '+' + key;
          break;
        case -1:
          ones += '-' + key;
          break;
        case 2:
          twos += '+' + key;
          break;
        case -2:
          twos += '-' + key;
          break;
        default:
          throw new Error('Inconsistent squaring table');
      }
    }
    let line = `res[${i}]=`;
    if (twos.length) {
      if (ones.length) {
        line += `${ones}+2*(${twos})`;
      } else {
        line += `2*(${twos})`;
      }
      squareLines.push(line);
    } else if (ones.length) {
      squareLines.push(line + ones);
    }
  }

  type binaryOp = (other: AlgebraElement) => AlgebraElement;
  const prelude = 'const res=new this.constructor();\nconst t=this;\n';
  const finale = '\nreturn res;';
  Result.prototype.add = new Function(
    'o',
    prelude + addInner + finale
  ) as binaryOp;
  Result.prototype.sub = new Function(
    'o',
    prelude + subInner + finale
  ) as binaryOp;
  Result.prototype.mul = new Function(
    'o',
    prelude + mulLines.join('\n') + finale
  ) as binaryOp;
  Result.prototype.wedge = new Function(
    'o',
    prelude + wedgeLines.join('\n') + finale
  ) as binaryOp;
  Result.prototype.vee = new Function(
    'o',
    prelude + veeLines.join('\n') + finale
  ) as binaryOp;
  Result.prototype.dot = new Function(
    'o',
    prelude + dotLines.join('\n') + finale
  ) as binaryOp;
  Result.prototype.dotL = new Function(
    'o',
    prelude + dotLeftLines.join('\n') + finale
  ) as binaryOp;

  Result.prototype.square = new Function(
    '',
    prelude + squareLines.join('\n') + finale
  ) as () => AlgebraElement;

  // We lose the option to negotiate numeric precision but gain speed
  Result.prototype.lmul = function (other: AlgebraElement) {
    return other.mul(this);
  };
  Result.prototype.lwedge = function (other: AlgebraElement) {
    return other.wedge(this);
  };
  Result.prototype.lvee = function (other: AlgebraElement) {
    return other.vee(this);
  };
  Result.prototype.ldotL = function (other: AlgebraElement) {
    return other.dotL(this);
  };

  return Result;
}

export default Algebra;
