// Float32Array-like
export declare class ElementBaseType {
  constructor(values?: number | Iterable<number>);
  [index: number]: number;
  get length(): number;
  fill(fillValue: number): ElementBaseType;
  [Symbol.iterator](): Iterator<number>;
}

export declare class AlgebraElement extends ElementBaseType {
  fill(fillValue: number): AlgebraElement;

  // Comparisons
  equals(other: AlgebraElement): boolean;
  closeTo(other: AlgebraElement): boolean;

  // Validation
  hasNaN(): boolean;
  hasInfinity(): boolean;

  // Getters / setters
  get s(): number; // Scalar part
  set s(value: number);
  get ps(): number; // Pseudoscalar part
  set ps(value: number);

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
  normalize(): AlgebraElement;

  // Scalar operations
  scale(scalar: number): AlgebraElement;

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
export function closeTo(a: AlgebraElement, b: AlgebraElement): boolean {
  return a.closeTo(b);
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
export function normalize(element: AlgebraElement): AlgebraElement {
  return element.normalize();
}

// Scalar operations
export function scale(element: AlgebraElement, scalar: number): AlgebraElement {
  return element.scale(scalar);
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

const MAX_DIMENSIONS = 36;

export default function Algebra(
  p: number,
  q = 0,
  r = 0,
  baseType: typeof ElementBaseType = Float32Array
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

  // Geometric product between basis vectors of index a and b
  function basisMul(a: number, b: number) {
    const aIndices = [];
    const bIndices = [];
    for (let i = 0; i < dimensions; ++i) {
      const p = 1 << i;
      if (p & a) {
        aIndices.push(i);
      }
      if (p & b) {
        bIndices.push(i);
      }
    }
    const indices = aIndices.concat(bIndices);
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

    fill(fillValue: number) {
      super.fill(fillValue);
      return this;
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
        result[i] = this[i] * [1, 1, -1, -1][bitCount(i) % 4];
      }
      return result;
    }

    involute(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i] * [1, -1, 1, -1][bitCount(i) % 4];
      }
      return result;
    }

    conjugate(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        result[i] = this[i] * [1, -1, -1, 1][bitCount(i) % 4];
      }
      return result;
    }

    // For all Ex = AlgebraClass.basisVector(...x)
    // Ex.mul(Ex.dual()) === AlgebraClass.pseudoscalar()
    dual(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] = this[i] * basisMul(i, dualIndex);
      }
      return result;
    }

    undual(): AlgebraElement {
      const result = new AlgebraClass();
      for (let i = 0; i < this.length; ++i) {
        const dualIndex = indexMask ^ i;
        result[dualIndex] = this[i] * basisMul(dualIndex, i);
      }
      return result;
    }

    normalize(): AlgebraElement {
      return this.scale(1 / this.norm());
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
          const normSquared = this.mul(this.conjugate());
          const n34 = normSquared.negateGrades(3, 4);
          return this.conjugate()
            .mul(n34)
            .scale(1 / normSquared.mul(n34).s);
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
          result[i ^ j] += this[i] * other[j] * basisMul(i, j);
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
          result[i ^ j] += this[i] * other[j] * basisMul(j, i);
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
            result[i ^ j] += this[i] * other[j] * basisMul(i, j);
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
            result[i ^ j] += this[i] * other[j] * basisMul(j, i);
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
            result[target] += this[i] * other[j] * basisMul(i, j);
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

    even() {
      const result = AlgebraClass.zero();
      for (let i = 0; i < this.length; ++i) {
        if (bitCount(i) % 2 === 0) {
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
      result[indices.map(i => 1 << i).reduce((a, b) => a | b, 0)] =
        sortSign(indices);
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

  return AlgebraClass;
}
