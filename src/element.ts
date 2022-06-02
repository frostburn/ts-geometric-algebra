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
  // See star() for forward implementation
  unstar(): AlgebraElement;
  hodge(): AlgebraElement;
  unhodge(): AlgebraElement;

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
  rmul(other: AlgebraElement): AlgebraElement;
  div(other: AlgebraElement): AlgebraElement;
  ldiv(other: AlgebraElement): AlgebraElement;
  ldivs(other: AlgebraElement): AlgebraElement;
  wedge(other: AlgebraElement): AlgebraElement;
  rwedge(other: AlgebraElement): AlgebraElement;
  vee(other: AlgebraElement): AlgebraElement;
  rvee(other: AlgebraElement): AlgebraElement;
  rotorMean(other: AlgebraElement): AlgebraElement;
  // Contractions
  contract(
    other: AlgebraElement,
    criterion: (r: number, s: number) => number
  ): AlgebraElement;
  dot(other: AlgebraElement): AlgebraElement; // Symmetric contraction
  dotL(other: AlgebraElement): AlgebraElement; // Left contraction
  dotR(other: AlgebraElement): AlgebraElement; // Right contraction
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
  plus(scalar: number): AlgebraElement;
  rescale(scalar: number): this;
  accumulate(other: AlgebraElement): this;
  split(iter?: number): AlgebraElement[];
  factorize(iter?: number): AlgebraElement[];

  // Construction
  static zero(): AlgebraElement;
  static scalar(magnitude?: number): AlgebraElement;
  static pseudoscalar(magnitude?: number): AlgebraElement;
  static basisVector(...indices: number[]): AlgebraElement;
  static fromVector(values: Iterable<number>, grade?: number): AlgebraElement;
  static fromRotor(values: Iterable<number>): AlgebraElement;
  static fromGanja(values: Iterable<number>): AlgebraElement;

  // Algebra information
  static get dimensions(): number;
  static get size(): number;
}
