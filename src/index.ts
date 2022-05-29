import {eigenValues} from './utils';

// Float32Array-like
export declare class ElementBaseType {
  constructor(values?: number | Iterable<number>);
  [index: number]: number;
  [Symbol.iterator](): Iterator<number>;

  copyWithin(target: number, start: number, end?: number): this;
  every(
    predicate: (
      value: number,
      index: number,
      array: ElementBaseType
    ) => unknown,
    thisArg?: any
  ): boolean;
  fill(value: number, start?: number, end?: number): this;
  filter(
    predicate: (value: number, index: number, array: ElementBaseType) => any,
    thisArg?: any
  ): ElementBaseType;
  find(
    predicate: (value: number, index: number, obj: ElementBaseType) => boolean,
    thisArg?: any
  ): number | undefined;
  findIndex(
    predicate: (value: number, index: number, obj: ElementBaseType) => boolean,
    thisArg?: any
  ): number;
  forEach(
    callbackfn: (value: number, index: number, array: ElementBaseType) => void,
    thisArg?: any
  ): void;
  indexOf(searchElement: number, fromIndex?: number): number;
  join(separator?: string): string;
  lastIndexOf(searchElement: number, fromIndex?: number): number;
  readonly length: number;
  map(
    callbackfn: (
      value: number,
      index: number,
      array: ElementBaseType
    ) => number,
    thisArg?: any
  ): ElementBaseType;
  reduce(
    callbackfn: (
      previousValue: number,
      currentValue: number,
      currentIndex: number,
      array: ElementBaseType
    ) => number
  ): number;
  reduce(
    callbackfn: (
      previousValue: number,
      currentValue: number,
      currentIndex: number,
      array: ElementBaseType
    ) => number,
    initialValue: number
  ): number;
  reduce<U>(
    callbackfn: (
      previousValue: U,
      currentValue: number,
      currentIndex: number,
      array: ElementBaseType
    ) => U,
    initialValue: U
  ): U;
  reduceRight(
    callbackfn: (
      previousValue: number,
      currentValue: number,
      currentIndex: number,
      array: ElementBaseType
    ) => number
  ): number;
  reduceRight(
    callbackfn: (
      previousValue: number,
      currentValue: number,
      currentIndex: number,
      array: ElementBaseType
    ) => number,
    initialValue: number
  ): number;
  reduceRight<U>(
    callbackfn: (
      previousValue: U,
      currentValue: number,
      currentIndex: number,
      array: ElementBaseType
    ) => U,
    initialValue: U
  ): U;
  reverse(): ElementBaseType;
  set(array: ArrayLike<number>, offset?: number): void;
  slice(start?: number, end?: number): ElementBaseType;
  some(
    predicate: (
      value: number,
      index: number,
      array: ElementBaseType
    ) => unknown,
    thisArg?: any
  ): boolean;
  sort(compareFn?: (a: number, b: number) => number): this;
  subarray(begin?: number, end?: number): ElementBaseType;
  toLocaleString(): string;
  toString(): string;
  valueOf(): ElementBaseType;
}

export declare class AlgebraElement extends ElementBaseType {
  // Comparisons
  equals(other: AlgebraElement): boolean;
  closeTo(other: AlgebraElement, tolerance?: number): boolean;

  // Validation
  hasNaN(): boolean;
  hasInfinity(): boolean;
  isNil(tolerance?: number): boolean;
  isGrade(grade: number, tolerance?: number): boolean;

  // Getters / setters
  get s(): number; // Scalar part
  set s(value: number);
  get ps(): number; // Pseudoscalar part
  set ps(value: number);

  getAt(...indices: number[]): number;
  setAt(...indicesAndValue: number[]): this;

  // Unary scalar operations
  norm(): number;
  vnorm(): number;

  // Unary operations
  neg(): AlgebraElement;
  cwAbs(): AlgebraElement;
  involute(): AlgebraElement;
  rev(): AlgebraElement;
  conjugate(): AlgebraElement;
  dual(): AlgebraElement;
  undual(): AlgebraElement;
  inverse(): AlgebraElement;
  normalize(newNorm?: number): AlgebraElement;
  exp(forceTaylor?: boolean, numTaylorTerms?: number): AlgebraElement;
  log(): AlgebraElement;
  clone(): AlgebraElement;

  // Scalar operations
  scale(scalar: number): AlgebraElement;
  pow(scalar: number): AlgebraElement;

  // Multi-scalar operations
  applyWeights(weights: number[]): AlgebraElement;

  // Index operations
  negateGrades(...grades: number[]): AlgebraElement;

  // Binary operations
  add(other: AlgebraElement): AlgebraElement;
  sub(other: AlgebraElement): AlgebraElement;
  mul(other: AlgebraElement): AlgebraElement;
  rmul(other: AlgebraElement): AlgebraElement;
  div(other: AlgebraElement): AlgebraElement;
  ldiv(other: AlgebraElement): AlgebraElement;
  ldivs(other: AlgebraElement): AlgebraElement;
  wedge(other: AlgebraElement): AlgebraElement;
  rwedge(other: AlgebraElement): AlgebraElement;
  vee(other: AlgebraElement): AlgebraElement;
  rvee(other: AlgebraElement): AlgebraElement;
  // Contractions
  contract(
    other: AlgebraElement,
    criterion: (r: number, s: number) => number
  ): AlgebraElement;
  dot(other: AlgebraElement): AlgebraElement; // Symmetric contraction
  dotL(other: AlgebraElement): AlgebraElement; // Left contraction
  dotR(other: AlgebraElement): AlgebraElement; // Right contraction
  star(other: AlgebraElement): AlgebraElement; // Scalar product

  // Subsets
  even(): AlgebraElement;
  grade(grade: number): AlgebraElement;

  // Deconstruction
  vector(grade?: number): ElementBaseType;
  ganja(): ElementBaseType;

  // Misc
  accumulate(other: AlgebraElement): this;
  split(iter?: number): AlgebraElement[];
  factorize(iter?: number): AlgebraElement[];

  // Construction
  static zero(): AlgebraElement;
  static scalar(magnitude?: number): AlgebraElement;
  static pseudoscalar(magnitude?: number): AlgebraElement;
  static basisVector(...indices: number[]): AlgebraElement;
  static fromVector(values: Iterable<number>, grade?: number): AlgebraElement;
  static fromGanja(values: Iterable<number>): AlgebraElement;

  // Algebra information
  static get dimensions(): number;
  static get size(): number;
}

// Comparisons using two arguments
export function equals(a: AlgebraElement, b: AlgebraElement): boolean {
  return a.equals(b);
}
export function closeTo(
  a: AlgebraElement,
  b: AlgebraElement,
  tolerance?: number
): boolean {
  return a.closeTo(b, tolerance);
}

// Validation
export function hasNaN(element: AlgebraElement): boolean {
  return element.hasNaN();
}
export function hasInfinity(element: AlgebraElement): boolean {
  return element.hasInfinity();
}
export function isNil(element: AlgebraElement, tolerance?: number): boolean {
  return element.isNil(tolerance);
}
export function isGrade(
  element: AlgebraElement,
  grade: number,
  tolerance?: number
): boolean {
  return element.isGrade(grade, tolerance);
}

// Unary scalar operations using one argument
export function norm(element: AlgebraElement): number {
  return element.norm();
}
export function vnorm(element: AlgebraElement): number {
  return element.vnorm();
}

// Unary operations using one argument
export function neg(element: AlgebraElement): AlgebraElement {
  return element.neg();
}
export function cwAbs(element: AlgebraElement): AlgebraElement {
  return element.cwAbs();
}
export function involute(element: AlgebraElement): AlgebraElement {
  return element.involute();
}
export function rev(element: AlgebraElement): AlgebraElement {
  return element.rev();
}
export function conjugate(element: AlgebraElement): AlgebraElement {
  return element.conjugate();
}
export function dual(element: AlgebraElement): AlgebraElement {
  return element.dual();
}
export function undual(element: AlgebraElement): AlgebraElement {
  return element.undual();
}
export function inverse(element: AlgebraElement): AlgebraElement {
  return element.inverse();
}
export function normalize(
  element: AlgebraElement,
  newNorm?: number
): AlgebraElement {
  return element.normalize(newNorm);
}
export function exp(
  element: AlgebraElement,
  forceTaylor?: boolean,
  numTaylorTerms?: number
): AlgebraElement {
  return element.exp(forceTaylor, numTaylorTerms);
}
export function clone(element: AlgebraElement) {
  return element.clone();
}

// Scalar operations
export function scale(element: AlgebraElement, scalar: number): AlgebraElement {
  return element.scale(scalar);
}
export function pow(element: AlgebraElement, power: number): AlgebraElement {
  return element.pow(power);
}

// Multi-scalar operations
export function applyWeights(
  element: AlgebraElement,
  weights: number[]
): AlgebraElement {
  return element.applyWeights(weights);
}

// Index operations
export function negateGrades(
  element: AlgebraElement,
  ...grades: number[]
): AlgebraElement {
  return element.negateGrades(...grades);
}

// Binary operations using two arguments
export function add(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.add(b);
}
export function sub(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.sub(b);
}
export function mul(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.mul(b);
}
export function div(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.div(b);
}
export function ldivs(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.ldivs(b);
}
export function wedge(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.wedge(b);
}
export function vee(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.vee(b);
}
// Contractions
export function contract(
  a: AlgebraElement,
  b: AlgebraElement,
  criterion: (r: number, s: number) => number
): AlgebraElement {
  return a.contract(b, criterion);
}
export function dot(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.dot(b);
}
export function dotL(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.dotL(b);
}
export function dotR(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.dotR(b);
}
export function star(a: AlgebraElement, b: AlgebraElement): AlgebraElement {
  return a.star(b);
}

// Subsets
export function even(element: AlgebraElement): AlgebraElement {
  return element.even();
}
export function grade(element: AlgebraElement, grade: number): AlgebraElement {
  return element.grade(grade);
}

// https://stackoverflow.com/a/43122214
function bitCount(n: number) {
  n = n - ((n >> 1) & 0x55555555);
  n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
  return (((n + (n >> 4)) & 0xf0f0f0f) * 0x1010101) >> 24;
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

export default function Algebra(
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
        if (Math.abs(this[i] - other[i]) > tolerance) {
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
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = -this[i];
      }
      return result;
    }

    cwAbs(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = Math.abs(this[i]);
      }
      return result;
    }

    rev(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = bitCount(i) & 2 ? -this[i] : this[i];
      }
      return result;
    }

    involute(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = bitCount(i) & 1 ? -this[i] : this[i];
      }
      return result;
    }

    conjugate(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = (bitCount(i) + 1) & 2 ? -this[i] : this[i];
      }
      return result;
    }

    // For all Ex = AlgebraClass.basisVector(...x)
    // Ex.mul(Ex.dual()) === AlgebraClass.pseudoscalar()
    dual(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] = this[i] * mulTable[i][dualIndex];
      }
      return result;
    }

    undual(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] = this[i] * mulTable[dualIndex][i];
      }
      return result;
    }

    normalize(newNorm = 1): AlgebraElement {
      return this.scale(newNorm / this.norm());
    }

    exp(forceTaylor = false, numTaylorTerms = 32) {
      if (!forceTaylor) {
        if (dimensions === 0) {
          return AlgebraClass.scalar(Math.exp(this.s));
        } else if (dimensions === 1) {
          const expS = Math.exp(this.s);
          if (p) {
            return new AlgebraClass([
              expS * Math.cosh(this.ps),
              expS * Math.sinh(this.ps),
            ]);
          } else if (q) {
            return new AlgebraClass([
              expS * Math.cos(this.ps),
              expS * Math.sin(this.ps),
            ]);
          } else if (r) {
            return new AlgebraClass([expS, expS * this.ps]);
          }
        } else if (dimensions === 2) {
          if (q === 2) {
            const expS = Math.exp(this.s);
            const imag = this.clone();
            imag.s = 0;
            const imagNorm = imag.vnorm();
            const result = imag.scale((expS * Math.sin(imagNorm)) / imagNorm);
            result.s = expS * Math.cos(imagNorm);
            return result;
          }
        }
      }
      if (!forceTaylor) {
        // Closed form exp
        const grade2 = this.clone();
        grade2.s = 0;
        if (grade2.isGrade(2)) {
          return grade2.split().reduce((total, simple) => {
            const square = simple.mul(simple).s,
              len = Math.sqrt(Math.abs(square));
            if (len <= 1e-5) {
              simple.s += 1;
            } else if (square < 0) {
              simple = simple.scale(Math.sin(len) / len);
              simple.s += Math.cos(len);
            } else {
              simple = simple.scale(Math.sinh(len) / len);
              simple.s += Math.cosh(len);
            }
            return total.mul(simple);
          }, AlgebraClass.scalar(Math.exp(this.s)));
        }
      }

      // Taylor series
      let result = AlgebraClass.scalar();
      let term = AlgebraClass.scalar();
      for (let i = 1; i < numTaylorTerms; ++i) {
        term = term.mul(this.scale(1 / i));
        result = result.add(term);
      }
      return result;
    }

    log() {
      if (dimensions === 0) {
        return AlgebraClass.scalar(Math.log(this.s));
      } else if (dimensions === 1) {
        if (p) {
          const norm = Math.sqrt(this.s ** 2 - this.ps ** 2);
          return new AlgebraClass([Math.log(norm), Math.asinh(this.ps / norm)]);
        } else if (q) {
          const norm = Math.hypot(this.s, this.ps);
          return new AlgebraClass([
            Math.log(norm),
            Math.atan2(this.ps, this.s),
          ]);
        } else if (r) {
          return new AlgebraClass([Math.log(this.s), this.ps / this.s]);
        }
      } else if (dimensions === 2) {
        if (q === 2) {
          const norm = this.vnorm();
          const imag = this.clone();
          imag.s = 0;
          const imagNorm = imag.vnorm();
          const result = imag.scale(Math.acos(this.s / norm) / imagNorm);
          result.s = Math.log(norm);
          return result;
        }
      }

      return this.factorize().reduce((sum, bi) => {
        const [ci, si] = [bi.s, bi.grade(2)];
        const square = si.mul(si).s;
        const len = Math.sqrt(Math.abs(square));
        if (Math.abs(square) < 1e-5) return sum.add(si);
        if (square < 0) return sum.add(si.scale(Math.acos(ci) / len));
        return sum.add(si.scale(Math.acosh(ci) / len));
      }, AlgebraClass.zero());
    }

    clone(): AlgebraElement {
      return new AlgebraClass(this);
    }

    negateGrades(...grades: number[]): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = grades.includes(bitCount(i)) ? -this[i] : this[i];
      }
      return result;
    }

    inverse(): AlgebraElement {
      // Matrix-free inverses up to 5D
      // http://repository.essex.ac.uk/17282/1/TechReport_CES-534.pdf
      switch (dimensions) {
        case 0:
          return AlgebraClass.scalar(1 / this.s);
        case 1:
          const involute = this.involute();
          return involute.scale(1 / this.mul(involute).s);
        case 2:
          const conjugate = this.conjugate();
          return conjugate.scale(1 / this.mul(conjugate).s);
        case 3:
          const reverse = this.rev();
          const involute3 = this.involute();
          const conjugate3 = this.conjugate();
          return reverse
            .mul(involute3)
            .mul(conjugate3)
            .scale(1 / this.mul(conjugate3).mul(involute3).mul(reverse).s);
        case 4:
          const modulus = this.mul(this.conjugate());
          const n34 = modulus.negateGrades(3, 4);
          return this.conjugate()
            .mul(n34)
            .scale(1 / modulus.mul(n34).s);
        case 5:
          const civ = this.conjugate().mul(this.involute()).mul(this.rev());
          const tciv = this.mul(civ);
          const tciv14 = tciv.negateGrades(1, 4);
          return civ.mul(tciv14).scale(1 / tciv.mul(tciv14).s);
        default:
          // Shirokov inverse
          const N = 1 << (((dimensions + 1) / 2) | 0);
          let Uk = this.scale(1);
          let adjU: AlgebraElement;
          for (let k = 1; k < N; ++k) {
            adjU = Uk.sub(AlgebraClass.scalar((N / k) * Uk.s));
            Uk = this.mul(adjU);
          }
          return Uk.s === 0 ? AlgebraClass.zero() : adjU!.scale(1 / Uk.s);
      }
    }

    scale(scalar: number): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i] * scalar;
      }
      return result;
    }

    pow(power: number): AlgebraElement {
      if (power !== Math.round(power)) {
        return this.log().scale(power).exp();
      }
      if (power === 0) {
        return AlgebraClass.scalar();
      }
      if (power === 1) {
        return this.clone();
      }
      if (power === 2) {
        return this.mul(this);
      }
      if (power > 0) {
        let result = AlgebraClass.scalar();
        let powerOfTwo = this.clone();
        while (power) {
          if (power & 1) {
            result = result.mul(powerOfTwo);
          }
          powerOfTwo = powerOfTwo.mul(powerOfTwo);
          power >>= 1;
        }
        return result;
      }
      return this.inverse().pow(-power);
    }

    applyWeights(weights: number[]) {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i];
        for (let j = 0; j < weights.length; ++j) {
          if (i & (1 << j)) {
            result[i] *= weights[j];
          }
        }
      }
      return result;
    }

    add(other: AlgebraElement): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i] + other[i];
      }
      return result;
    }

    sub(other: AlgebraElement): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i] - other[i];
      }
      return result;
    }

    mul(other: AlgebraElement): AlgebraElement {
      const result = AlgebraClass.zero();
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

    rmul(other: AlgebraElement): AlgebraElement {
      const result = AlgebraClass.zero();
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
      const result = AlgebraClass.zero();
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

    rwedge(other: AlgebraElement): AlgebraElement {
      const result = AlgebraClass.zero();
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

    rvee(other: AlgebraElement): AlgebraElement {
      return other.dual().rwedge(this.dual()).undual();
    }

    contract(
      other: AlgebraElement,
      criterion: (r: number, s: number) => number
    ): AlgebraElement {
      const result = AlgebraClass.zero();
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

    dotR(other: AlgebraElement) {
      return this.contract(other, right);
    }

    star(other: AlgebraElement) {
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

    even() {
      const result = AlgebraClass.zero();
      for (let i = 0; i < this.length; ++i) {
        if (!(bitCount(i) & 1)) {
          result[i] = this[i];
        }
      }
      return result;
    }

    grade(grade: number) {
      const result = AlgebraClass.zero();
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

    // Ganja.js compatible representation
    ganja() {
      return new baseType(indexString.map(g => this[g[0]]));
    }

    accumulate(other: AlgebraElement) {
      for (let i = 0; i < this.length; ++i) {
        this[i] += other[i];
      }
      return this;
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
        const D = 0.5 * Math.sqrt(TDT * TDT - TWT.mul(TWT).s);
        eigen = [0.5 * TDT + D, 0.5 * TDT - D].sort(
          (a, b) => Math.abs(a) - Math.abs(b)
        );
      } else {
        // For >6D, closed form solutions of the characteristic polyn. are impossible, use eigenvalues of companion matrix.
        const Wis = Wi.map((W, i) => W.mul(W).s * (-1) ** (k - i + (k % 2)));
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
      Wi = [AlgebraClass.scalar(), ...Wi, AlgebraClass.zero()];
      const sum = AlgebraClass.zero();
      const k2 = Math.floor(k / 2);
      const res: AlgebraElement[] = eigen.slice(1).map(v => {
        const N = AlgebraClass.zero();
        const DN = AlgebraClass.zero();
        for (let i = 0; i <= k2; ++i) {
          N.accumulate(Wi[2 * i + 1].scale(v ** (k2 - i)));
          DN.accumulate(Wi[2 * i].scale(v ** (k2 - i)));
        }
        if (DN.vnorm() === 0) return AlgebraClass.zero();
        const ret = N.div(DN);
        sum.accumulate(ret);
        return ret;
      });
      return [this.sub(sum), ...res]; // Smallest eigvalue becomes B-rest
    }

    // Factorize a motor
    factorize(iter = 50) {
      const S = this.grade(2).split(iter);
      const R = S.slice(0, S.length - 1).map(Mi => {
        Mi.s += this.s;
        const scale = Math.sqrt(Mi.rev().mul(Mi).s);
        return Mi.scale(1 / scale);
      });
      R.push(
        R.reduce((tot, fact) => tot.mul(fact.rev()), AlgebraClass.scalar()).mul(
          this
        )
      );
      return R;
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

    static basisVector(...indices: number[]): AlgebraElement {
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

  if (!unroll) {
    return AlgebraClass;
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

  type binaryOp = (other: AlgebraElement) => AlgebraElement;
  const prelude = 'const res=new this.constructor();\nconst t=this;\n';
  const finale = '\nreturn res;';
  AlgebraClass.prototype.add = new Function(
    'o',
    prelude + addInner + finale
  ) as binaryOp;
  AlgebraClass.prototype.sub = new Function(
    'o',
    prelude + subInner + finale
  ) as binaryOp;
  AlgebraClass.prototype.mul = new Function(
    'o',
    prelude + mulLines.join('\n') + finale
  ) as binaryOp;
  AlgebraClass.prototype.wedge = new Function(
    'o',
    prelude + wedgeLines.join('\n') + finale
  ) as binaryOp;
  AlgebraClass.prototype.vee = new Function(
    'o',
    prelude + veeLines.join('\n') + finale
  ) as binaryOp;
  AlgebraClass.prototype.dot = new Function(
    'o',
    prelude + dotLines.join('\n') + finale
  ) as binaryOp;
  AlgebraClass.prototype.dotL = new Function(
    'o',
    prelude + dotLeftLines.join('\n') + finale
  ) as binaryOp;

  // We lose the option to negotiate numeric precision but gain speed
  AlgebraClass.prototype.rmul = function (other: AlgebraElement) {
    return other.mul(this);
  };
  AlgebraClass.prototype.rwedge = function (other: AlgebraElement) {
    return other.wedge(this);
  };
  AlgebraClass.prototype.rvee = function (other: AlgebraElement) {
    return other.vee(this);
  };

  return AlgebraClass;
}
