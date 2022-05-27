export type ElementBaseType = typeof Float32Array | typeof Float64Array;

export declare class AlgebraElement {
  constructor(values?: Iterable<number>);
  [index: number]: number;

  // Float32Array methods
  get length(): number;
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
  dot(other: AlgebraElement): AlgebraElement;
  wedge(other: AlgebraElement): AlgebraElement;
  rwedge(other: AlgebraElement): AlgebraElement;
  vee(other: AlgebraElement): AlgebraElement;
  rvee(other: AlgebraElement): AlgebraElement;

  // Subsets
  even(): AlgebraElement;
  grade(grade: number): AlgebraElement;

  // Deconstruction
  vector(grade?: number): Float32Array | Float64Array;
  ganja(): Float32Array | Float64Array;

  // Construction
  static zero(): AlgebraElement;
  static scalar(magnitude?: number): AlgebraElement;
  static pseudoscalar(magnitude?: number): AlgebraElement;
  static basisVector(...indices: number[]): AlgebraElement;
  static fromVector(values: Iterable<number>, grade?: number): AlgebraElement;
  static fromGanja(values: Iterable<number>): AlgebraElement;

  // Binary operations using two arguments
  static add(a: AlgebraElement, b: AlgebraElement): AlgebraElement;
  static sub(a: AlgebraElement, b: AlgebraElement): AlgebraElement;
  static mul(a: AlgebraElement, b: AlgebraElement): AlgebraElement;
  static div(a: AlgebraElement, b: AlgebraElement): AlgebraElement;
  static dot(a: AlgebraElement, b: AlgebraElement): AlgebraElement;
  static wedge(a: AlgebraElement, b: AlgebraElement): AlgebraElement;
  static vee(a: AlgebraElement, b: AlgebraElement): AlgebraElement;

  // Algebra information
  static get dimensions(): number;
  static get size(): number;
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

export default function Algebra(
  p: number,
  q = 0,
  r = 0,
  baseType: ElementBaseType = Float32Array
): typeof AlgebraElement {
  const metric: number[] = [];
  for (let i = 0; i < p; ++i) {
    metric.push(1);
  }
  for (let i = 0; i < q; ++i) {
    metric.push(-1);
  }
  for (let i = 0; i < r; ++i) {
    metric.push(0);
  }
  const dimensions = p + q + r;
  const size = 1 << dimensions;
  const indexMask = size - 1;

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
        str += j.toString();
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

  // XXX: That `baseType as typeof Float32Array` shouldn't be necessary,
  // but VSCode keeps complaining.
  class AlgebraClass extends (baseType as typeof Float32Array) {
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
          return this.involute().scale(1 / this.mul(this.involute()).s);
        case 2:
          return this.conjugate().scale(1 / this.mul(this.conjugate()).s);
        case 3:
          return this.rev()
            .mul(this.involute())
            .mul(this.conjugate())
            .scale(
              1 /
                this.mul(this.conjugate()).mul(this.involute()).mul(this.rev())
                  .s
            );
        case 4:
          return this.conjugate()
            .mul(this.mul(this.conjugate()).negateGrades(3, 4))
            .scale(
              1 /
                this.mul(this.conjugate()).mul(
                  this.mul(this.conjugate()).negateGrades(3, 4)
                ).s
            );
        case 5:
          return this.conjugate()
            .mul(this.involute())
            .mul(this.rev())
            .mul(
              this.mul(this.conjugate())
                .mul(this.involute())
                .mul(this.rev())
                .negateGrades(1, 4)
            )
            .scale(
              1 /
                this.mul(this.conjugate())
                  .mul(this.involute())
                  .mul(this.rev())
                  .mul(
                    this.mul(this.conjugate())
                      .mul(this.involute())
                      .mul(this.rev())
                      .negateGrades(1, 4)
                  ).s
            );
        default:
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

    dot(other: AlgebraElement): AlgebraElement {
      const result = AlgebraClass.zero();
      for (let i = 0; i < this.length; ++i) {
        for (let j = 0; j < other.length; ++j) {
          result[i ^ j] +=
            this[i] * other[j] * (mulTable[i][j] + mulTable[j][i]) * 0.5;
        }
      }
      return result;
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
      return this.dual().wedge(other.dual()).undual();
    }

    rvee(other: AlgebraElement): AlgebraElement {
      return this.dual().rwedge(other.dual()).undual();
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
        if (bitCount(i) === grade) {
          result.push(this[i]);
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
          if (bitCount(i) === grade) {
            result[i] = component;
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

    static add(a: AlgebraElement, b: AlgebraElement) {
      return a.add(b);
    }

    static sub(a: AlgebraElement, b: AlgebraElement) {
      return a.sub(b);
    }

    static mul(a: AlgebraElement, b: AlgebraElement) {
      return a.mul(b);
    }

    static div(a: AlgebraElement, b: AlgebraElement) {
      return a.div(b);
    }

    static dot(a: AlgebraElement, b: AlgebraElement) {
      return a.dot(b);
    }

    static wedge(a: AlgebraElement, b: AlgebraElement) {
      return a.wedge(b);
    }

    static vee(a: AlgebraElement, b: AlgebraElement) {
      return a.vee(b);
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
