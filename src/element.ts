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
  inverse(): AlgebraElement;
  square(): AlgebraElement;
  normalize(newNorm?: number): AlgebraElement;
  rotorNormalize(): AlgebraElement;
  sqrt(forceBabylon?: boolean, numIter?: number): AlgebraElement;
  rotorSqrt(): AlgebraElement;
  exp(forceTaylor?: boolean, numTaylorTerms?: number): AlgebraElement;
  bivectorExp(): AlgebraElement;
  log(): AlgebraElement;
  rotorLog(): AlgebraElement;
  clone(): AlgebraElement;
  // Dual Zoo
  dual(): AlgebraElement;
  undual(): AlgebraElement;
  podge(): AlgebraElement;
  unpodge(): AlgebraElement;
  podgeL(): AlgebraElement;
  unpodgeL(): AlgebraElement;
  // See star() for forward implementation
  unstar(): AlgebraElement;
  starL(): AlgebraElement;
  unstarL(): AlgebraElement;
  hodge(): AlgebraElement;
  unhodge(): AlgebraElement;
  hodgeL(): AlgebraElement;
  unhodgeL(): AlgebraElement;

  // Scalar operations
  scale(scalar: number): AlgebraElement;
  pow(scalar: number, splitStages?: number): AlgebraElement;

  // Multi-scalar operations
  applyWeights(weights: number[]): AlgebraElement;

  // Index operations
  negateGrades(...grades: number[]): AlgebraElement;

  // Binary operations
  add(other: AlgebraElement): AlgebraElement;
  sub(other: AlgebraElement): AlgebraElement;
  mul(other: AlgebraElement): AlgebraElement;
  lmul(other: AlgebraElement): AlgebraElement;
  div(other: AlgebraElement): AlgebraElement;
  ldiv(other: AlgebraElement): AlgebraElement;
  ldivs(other: AlgebraElement): AlgebraElement;
  wedge(other: AlgebraElement): AlgebraElement;
  lwedge(other: AlgebraElement): AlgebraElement;
  vee(other: AlgebraElement): AlgebraElement;
  lvee(other: AlgebraElement): AlgebraElement;
  rotorMean(other: AlgebraElement): AlgebraElement;
  // Contractions
  contract(
    other: AlgebraElement,
    criterion: (r: number, s: number) => number
  ): AlgebraElement;
  dot(other: AlgebraElement): AlgebraElement; // Symmetric contraction
  dotL(other: AlgebraElement): AlgebraElement; // Left contraction
  ldotL(other: AlgebraElement): AlgebraElement; // Left contraction
  dotR(other: AlgebraElement): AlgebraElement; // Right contraction
  ldotR(other: AlgebraElement): AlgebraElement; // Right contraction
  star(): AlgebraElement; // Dischord dual
  star(other: AlgebraElement): AlgebraElement; // Scalar product

  // Subsets
  imag(): AlgebraElement;
  even(): AlgebraElement;
  grade(grade: number): AlgebraElement;

  // Deconstruction
  vector(grade?: number): ElementBaseType;
  rotor(): ElementBaseType;
  ganja(): ElementBaseType;

  // Misc
  grades(threshold?: number): number[];
  plus(scalar: number): AlgebraElement;
  rescale(scalar: number): this;
  accumulate(other: AlgebraElement): this;
  split(iter?: number): AlgebraElement[];
  factorize(iter?: number): AlgebraElement[];

  // Construction
  static zero(): AlgebraElement;
  static scalar(magnitude?: number): AlgebraElement;
  static pseudoscalar(magnitude?: number): AlgebraElement;
  static basisBlade(...indices: number[]): AlgebraElement;
  static fromVector(values: Iterable<number>, grade?: number): AlgebraElement;
  static fromRotor(values: Iterable<number>): AlgebraElement;
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
export function inverse(element: AlgebraElement): AlgebraElement {
  return element.inverse();
}
export function square(element: AlgebraElement): AlgebraElement {
  return element.square();
}
export function normalize(
  element: AlgebraElement,
  newNorm?: number
): AlgebraElement {
  return element.normalize(newNorm);
}
export function rotorNormalize(element: AlgebraElement): AlgebraElement {
  return element.rotorNormalize();
}
export function sqrt(element: AlgebraElement): AlgebraElement {
  return element.sqrt();
}
export function rotorSqrt(element: AlgebraElement): AlgebraElement {
  return element.rotorSqrt();
}
export function exp(
  element: AlgebraElement,
  forceTaylor?: boolean,
  numTaylorTerms?: number
): AlgebraElement {
  return element.exp(forceTaylor, numTaylorTerms);
}
export function bivectorExp(element: AlgebraElement): AlgebraElement {
  return element.bivectorExp();
}
export function log(element: AlgebraElement): AlgebraElement {
  return element.log();
}
export function rotorLog(element: AlgebraElement): AlgebraElement {
  return element.rotorLog();
}
export function clone(element: AlgebraElement) {
  return element.clone();
}
// Dual Zoo
export function dual(element: AlgebraElement): AlgebraElement {
  return element.dual();
}
export function undual(element: AlgebraElement): AlgebraElement {
  return element.undual();
}
export function podge(element: AlgebraElement): AlgebraElement {
  return element.podge();
}
export function unpodge(element: AlgebraElement): AlgebraElement {
  return element.unpodge();
}
export function podgeL(element: AlgebraElement): AlgebraElement {
  return element.podgeL();
}
export function unpodgeL(element: AlgebraElement): AlgebraElement {
  return element.unpodgeL();
}
// See star overload for forward implementation
export function unstar(element: AlgebraElement): AlgebraElement {
  return element.unstar();
}
export function starL(element: AlgebraElement): AlgebraElement {
  return element.starL();
}
export function unstarL(element: AlgebraElement): AlgebraElement {
  return element.unstarL();
}
export function hodge(element: AlgebraElement): AlgebraElement {
  return element.hodge();
}
export function unhodge(element: AlgebraElement): AlgebraElement {
  return element.unhodge();
}
export function hodgeL(element: AlgebraElement): AlgebraElement {
  return element.hodgeL();
}
export function unhodgeL(element: AlgebraElement): AlgebraElement {
  return element.unhodgeL();
}

// Scalar operations
export function scale(element: AlgebraElement, scalar: number): AlgebraElement {
  return element.scale(scalar);
}
export function pow(
  element: AlgebraElement,
  power: number,
  splitStages?: number
): AlgebraElement {
  return element.pow(power, splitStages);
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
export function star(element: AlgebraElement): AlgebraElement;
export function star(a: AlgebraElement, b: AlgebraElement): AlgebraElement;
export function star(a: AlgebraElement, b?: AlgebraElement): AlgebraElement {
  if (b === undefined) {
    return a.star();
  }
  return a.star(b);
}

// Subsets
export function imag(element: AlgebraElement): AlgebraElement {
  return element.imag();
}
export function even(element: AlgebraElement): AlgebraElement {
  return element.even();
}
export function grade(element: AlgebraElement, grade: number): AlgebraElement {
  return element.grade(grade);
}

// Misc
export function grades(element: AlgebraElement, threshold?: number): number[] {
  return element.grades(threshold);
}
