// Specialized implementations for particular metrics

import {AlgebraElement} from './element';
import {
  Bivariate2D,
  complexSqrt,
  complexLog,
  complexExp,
  splitComplexSqrt,
  splitComplexExp,
  splitComplexLog,
  dualSqrt,
  dualExp,
  dualLog,
  sinc,
  sinch,
} from './utils';

function make2D(
  baseClass: typeof AlgebraElement,
  sqrt: Bivariate2D,
  exp: Bivariate2D,
  log: Bivariate2D,
  name: string
) {
  class Number2D extends baseClass {
    get algebra() {
      return Number2D;
    }

    sqrt(forceIter = false, numIter = 16): AlgebraElement {
      if (forceIter) {
        return super.sqrt(forceIter, numIter);
      }
      return new this.algebra(sqrt(this.s, this.ps));
    }

    exp(forceSeries = false, numTerms = 16): AlgebraElement {
      if (forceSeries) {
        return super.exp(forceSeries, numTerms);
      }
      return new this.algebra(exp(this.s, this.ps));
    }

    log(): AlgebraElement {
      return new this.algebra(log(this.s, this.ps));
    }

    inverse(): AlgebraElement {
      const involute = this.involute();
      return involute.scale(1 / this.mul(involute).s);
    }

    static zero() {
      return new Number2D().fill(0);
    }
    static scalar(magnitude = 1): AlgebraElement {
      return new Number2D(baseClass.scalar(magnitude));
    }
    static pseudoscalar(magnitude = 1): AlgebraElement {
      return new Number2D(baseClass.pseudoscalar(magnitude));
    }
    static basisBlade(...indices: number[]): AlgebraElement {
      return new Number2D(baseClass.basisBlade(...indices));
    }
    static fromVector(
      values: Iterable<number>,
      grade?: number
    ): AlgebraElement {
      return new Number2D(baseClass.fromVector(values, grade));
    }
    static fromRotor(values: Iterable<number>): AlgebraElement {
      return new Number2D(baseClass.fromRotor(values));
    }
    static fromGanja(values: Iterable<number>) {
      return new Number2D(baseClass.fromGanja(values));
    }
  }
  Object.defineProperty(Number2D, 'name', {value: name});
  return Number2D;
}

function makeSplitQuaternion(
  baseClass: typeof AlgebraElement,
  p1: number,
  p2: number,
  q: number
) {
  class SplitQuaternion extends baseClass {
    get algebra() {
      return SplitQuaternion;
    }

    imagNorm2() {
      return this[p1] ** 2 + this[p2] ** 2 - this[q] ** 2;
    }

    sqrt(forceIter = false, numIter = 16): AlgebraElement {
      if (forceIter) {
        return super.sqrt(forceIter, numIter);
      }
      const result = this.imag();
      const imagNorm2 = this.imagNorm2();
      const imagNorm = Math.sqrt(Math.abs(imagNorm2));
      let x, y;
      if (imagNorm2 < 0) {
        [x, y] = complexSqrt(this.s, imagNorm);
      } else {
        [x, y] = splitComplexSqrt(this.s, imagNorm);
      }
      if (imagNorm < 1e-5) {
        result.rescale(0.5);
      } else {
        result.rescale(y / imagNorm);
      }
      result.s = x;
      return result;
    }

    exp(forceSeries = false, numTerms = 16): AlgebraElement {
      if (forceSeries) {
        return super.exp(forceSeries, numTerms);
      }
      const expS = Math.exp(this.s);
      const result = this.imag();
      const imagNorm2 = this.imagNorm2();
      const imagNorm = Math.sqrt(Math.abs(imagNorm2));
      if (imagNorm < 1e-5) {
        result.rescale(expS);
      } else if (imagNorm2 < 0) {
        result.rescale((expS * Math.sin(imagNorm)) / imagNorm);
      } else {
        result.rescale((expS * Math.sinh(imagNorm)) / imagNorm);
      }
      if (imagNorm2 < 0) {
        result.s = expS * Math.cos(imagNorm);
      } else {
        result.s = expS * Math.cosh(imagNorm);
      }
      return result;
    }

    log(): AlgebraElement {
      const result = this.imag();
      const imagNorm2 = this.imagNorm2();
      const imagNorm = Math.sqrt(Math.abs(imagNorm2));
      const norm = Math.sqrt(this.s ** 2 - imagNorm2);
      if (imagNorm2 < 0) {
        result.rescale(Math.atan2(imagNorm, this.s) / imagNorm);
      } else {
        result.rescale(Math.asinh(imagNorm / norm) / imagNorm);
      }
      result.s = Math.log(norm);
      return result;
    }

    inverse(): AlgebraElement {
      const conjugate = this.conjugate();
      return conjugate.scale(1 / this.mul(conjugate).s);
    }

    static zero() {
      return new SplitQuaternion().fill(0);
    }
    static scalar(magnitude = 1): AlgebraElement {
      return new SplitQuaternion(baseClass.scalar(magnitude));
    }
    static pseudoscalar(magnitude = 1): AlgebraElement {
      return new SplitQuaternion(baseClass.pseudoscalar(magnitude));
    }
    static basisBlade(...indices: number[]): AlgebraElement {
      return new SplitQuaternion(baseClass.basisBlade(...indices));
    }
    static fromVector(
      values: Iterable<number>,
      grade?: number
    ): AlgebraElement {
      return new SplitQuaternion(baseClass.fromVector(values, grade));
    }
    static fromRotor(values: Iterable<number>): AlgebraElement {
      return new SplitQuaternion(baseClass.fromRotor(values));
    }
    static fromGanja(values: Iterable<number>) {
      return new SplitQuaternion(baseClass.fromGanja(values));
    }
  }
  return SplitQuaternion;
}

function makePGA1D(
  baseClass: typeof AlgebraElement,
  sqrt: Bivariate2D,
  co: (x: number) => number,
  sic: (x: number) => number,
  log: Bivariate2D,
  name: string,
  r = 1,
  pq = 2
) {
  class PGA1D extends baseClass {
    get algebra() {
      return PGA1D;
    }

    sqrt(forceIter = false, numIter = 16): AlgebraElement {
      if (forceIter) {
        return super.sqrt(forceIter, numIter);
      }
      const [x, y] = sqrt(this.s, this[pq]);
      const result = this.clone();
      result.s = x;
      result[r] *= 0.5 / x;
      result[pq] = y;
      result[3] *= 0.5 / x;
      return result;
    }

    exp(forceSeries = false, numTerms = 16): AlgebraElement {
      if (forceSeries) {
        return super.exp(forceSeries, numTerms);
      }
      const expS = Math.exp(this.s);
      const eSic = expS * sic(this[pq]);
      const result = this.clone();
      result.s = expS * co(this[pq]);
      result[r] *= eSic;
      result[pq] *= eSic;
      result[3] *= eSic;
      return result;
    }

    log(): AlgebraElement {
      const [x, y] = log(this.s, this[pq]);
      const result = this.clone();
      result.s = x;
      result[r] = (result[r] / this[pq]) * y;
      result[pq] = y;
      result[3] = (result[3] / this[pq]) * y;
      return result;
    }

    inverse(): AlgebraElement {
      const conjugate = this.conjugate();
      return conjugate.scale(1 / this.mul(conjugate).s);
    }

    static zero() {
      return new PGA1D().fill(0);
    }
    static scalar(magnitude = 1): AlgebraElement {
      return new PGA1D(baseClass.scalar(magnitude));
    }
    static pseudoscalar(magnitude = 1): AlgebraElement {
      return new PGA1D(baseClass.pseudoscalar(magnitude));
    }
    static basisBlade(...indices: number[]): AlgebraElement {
      return new PGA1D(baseClass.basisBlade(...indices));
    }
    static fromVector(
      values: Iterable<number>,
      grade?: number
    ): AlgebraElement {
      return new PGA1D(baseClass.fromVector(values, grade));
    }
    static fromRotor(values: Iterable<number>): AlgebraElement {
      return new PGA1D(baseClass.fromRotor(values));
    }
    static fromGanja(values: Iterable<number>) {
      return new PGA1D(baseClass.fromGanja(values));
    }
  }
  Object.defineProperty(PGA1D, 'name', {value: name});
  return PGA1D;
}

function makeSplitQuaternionPGA(
  baseClass: typeof AlgebraElement,
  p1: number,
  p2: number,
  q: number,
  sign: number
) {
  class SplitQuaternionPGA extends baseClass {
    get algebra() {
      return SplitQuaternionPGA;
    }

    imagNorm2() {
      return this[p1] * this[p1] + this[p2] * this[p2] - this[q] * this[q];
    }

    sqrt(forceIter = false, numIter = 16): AlgebraElement {
      if (forceIter) {
        return super.sqrt(forceIter, numIter);
      }

      const result = this.imag();
      const imagNorm2 = this.imagNorm2();
      const imagNorm = Math.sqrt(Math.abs(imagNorm2));
      let x, y;
      if (imagNorm2 < 0) {
        [x, y] = complexSqrt(this.s, imagNorm);
      } else {
        [x, y] = splitComplexSqrt(this.s, imagNorm);
      }
      result.s = x;
      let im = 0.5;
      if (imagNorm > 1e-5) {
        im = y / imagNorm;
      }
      result[2] *= im;
      result[4] *= im;
      result[6] *= im;

      const denom =
        x * x -
        result[p1] * result[p1] -
        result[p2] * result[p2] +
        result[q] * result[q];
      const ps =
        (-0.5 *
          (this[1] * result[6] +
            this[3] * result[4] -
            this[5] * result[2] -
            this.ps * x)) /
        denom;

      result[1] = (0.5 * this[1] + sign * result[6] * ps) / x;
      result[3] = (0.5 * this[3] - sign * result[4] * ps) / x;
      result[5] = (0.5 * this[5] + result[2] * ps) / x;
      result.ps = ps;

      return result;
    }

    inverse(): AlgebraElement {
      const reverse = this.rev();
      const involute = this.involute();
      const conjugate = this.conjugate();
      return reverse
        .mul(involute)
        .mul(conjugate)
        .scale(1 / this.mul(conjugate).mul(involute).mul(reverse).s);
    }

    static zero() {
      return new SplitQuaternionPGA().fill(0);
    }
    static scalar(magnitude = 1): AlgebraElement {
      return new SplitQuaternionPGA(baseClass.scalar(magnitude));
    }
    static pseudoscalar(magnitude = 1): AlgebraElement {
      return new SplitQuaternionPGA(baseClass.pseudoscalar(magnitude));
    }
    static basisBlade(...indices: number[]): AlgebraElement {
      return new SplitQuaternionPGA(baseClass.basisBlade(...indices));
    }
    static fromVector(
      values: Iterable<number>,
      grade?: number
    ): AlgebraElement {
      return new SplitQuaternionPGA(baseClass.fromVector(values, grade));
    }
    static fromRotor(values: Iterable<number>): AlgebraElement {
      return new SplitQuaternionPGA(baseClass.fromRotor(values));
    }
    static fromGanja(values: Iterable<number>) {
      return new SplitQuaternionPGA(baseClass.fromGanja(values));
    }
  }
  return SplitQuaternionPGA;
}

// https://arxiv.org/abs/2003.06873
function make3D(
  baseClass: typeof AlgebraElement,
  p1: number,
  p2: number,
  p3: number,
  q1: number,
  q2: number,
  q3: number,
  name: string
) {
  class Clifford3D extends baseClass {
    get algebra() {
      return Clifford3D;
    }

    sqrtScalars(signT = 1, signS = 1): [number, number] {
      const bS =
        this[0] * this[0] -
        this[p1] * this[p1] -
        this[p2] * this[p2] -
        this[p3] * this[p3] +
        this[q1] * this[q1] +
        this[q2] * this[q2] +
        this[q3] * this[q3] -
        this[7] * this[7];
      const bI =
        2 *
        (this[p3] * this[q1] -
          this[p2] * this[q2] +
          this[p1] * this[q3] -
          this[0] * this[7]);
      const sD = Math.hypot(bS, bI);
      let t, T;
      if (sD > bS) {
        t = 0.25 * (this[7] + signT * Math.sqrt(0.5 * (sD - bS)));
        T = 0.25 * ((signT * bI) / Math.sqrt(2 * (sD - bS)) - this[0]);
      } else {
        t = 0.25 * this[7];
        T = 0.25 * (signT * Math.sqrt(bS) - this[0]);
      }
      const t2 = Math.hypot(t, T);
      const s = signS * Math.sqrt(t2 - T);
      const S = t / s;
      return [s, S];
    }

    sqrt(forceIter = false, numIter = 16): AlgebraElement {
      if (forceIter) {
        return super.sqrt(forceIter, numIter);
      }
      const b1 = this[p1];
      const b2 = this[p2];
      const b3 = this[p3];
      const b12 = this[q1];
      const b13 = this[q2];
      const b23 = this[q3];
      const [s, S] = this.sqrtScalars();
      const s2 = 0.5 / (s * s + S * S);
      const v1 = (b1 * s + b23 * S) * s2;
      const v2 = (b2 * s - b13 * S) * s2;
      const v3 = (b3 * s + b12 * S) * s2;
      const V1 = (b23 * s - b1 * S) * s2;
      const V2 = -(b13 * s + b2 * S) * s2;
      const V3 = (b12 * s - b3 * S) * s2;

      const coords = Array(8).fill(0);
      coords[0] = s;
      coords[p1] = v1;
      coords[p2] = v2;
      coords[p3] = v3;
      coords[q1] = V3;
      coords[q2] = -V2;
      coords[q3] = V1;
      coords[7] = S;

      return new this.algebra(coords);
    }

    inverse(): AlgebraElement {
      const reverse = this.rev();
      const involute = this.involute();
      const conjugate = this.conjugate();
      return reverse
        .mul(involute)
        .mul(conjugate)
        .scale(1 / this.mul(conjugate).mul(involute).mul(reverse).s);
    }

    static zero() {
      return new Clifford3D().fill(0);
    }
    static scalar(magnitude = 1): AlgebraElement {
      return new Clifford3D(baseClass.scalar(magnitude));
    }
    static pseudoscalar(magnitude = 1): AlgebraElement {
      return new Clifford3D(baseClass.pseudoscalar(magnitude));
    }
    static basisBlade(...indices: number[]): AlgebraElement {
      return new Clifford3D(baseClass.basisBlade(...indices));
    }
    static fromVector(
      values: Iterable<number>,
      grade?: number
    ): AlgebraElement {
      return new Clifford3D(baseClass.fromVector(values, grade));
    }
    static fromRotor(values: Iterable<number>): AlgebraElement {
      return new Clifford3D(baseClass.fromRotor(values));
    }
    static fromGanja(values: Iterable<number>) {
      return new Clifford3D(baseClass.fromGanja(values));
    }
  }

  Object.defineProperty(Clifford3D, 'name', {value: name});

  return Clifford3D;
}

// https://www.researchgate.net/publication/360528787_Normalization_Square_Roots_and_the_Exponential_and_Logarithmic_Maps_in_Geometric_Algebras_of_Less_than_6D
export function pqrMixin(
  p: number,
  q: number,
  r: number,
  baseClass: typeof AlgebraElement
): typeof AlgebraElement {
  if (p === 0 && q === 0 && r === 0) {
    class Scalar extends baseClass {
      get algebra() {
        return Scalar;
      }

      sqrt(forceIter = false, numIter = 16): AlgebraElement {
        if (forceIter) {
          return super.sqrt(forceIter, numIter);
        }
        return this.algebra.scalar(Math.sqrt(this.s));
      }

      exp(forceSeries = false, numTerms = 16): AlgebraElement {
        if (forceSeries) {
          return super.exp(forceSeries, numTerms);
        }
        return this.algebra.scalar(Math.exp(this.s));
      }

      log(): AlgebraElement {
        return this.algebra.scalar(Math.log(this.s));
      }

      inverse(): AlgebraElement {
        return this.algebra.scalar(1 / this.s);
      }

      static zero() {
        return new Scalar().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Scalar(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Scalar(baseClass.pseudoscalar(magnitude));
      }
      static basisBlade(...indices: number[]): AlgebraElement {
        return new Scalar(baseClass.basisBlade(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Scalar(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Scalar(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Scalar(baseClass.fromGanja(values));
      }
    }
    return Scalar;
  }
  if (p === 1 && q === 0 && r === 0) {
    const SplitComplex = make2D(
      baseClass,
      splitComplexSqrt,
      splitComplexExp,
      splitComplexLog,
      'SplitComplex'
    );
    return SplitComplex;
  }
  if (p === 0 && q === 1 && r === 0) {
    const Complex = make2D(
      baseClass,
      complexSqrt,
      complexExp,
      complexLog,
      'Complex'
    );
    return Complex;
  }
  if (p === 0 && q === 0 && r === 1) {
    const Dual = make2D(baseClass, dualSqrt, dualExp, dualLog, 'Dual');
    return Dual;
  }
  if (p === 2 && q === 0 && r === 0) {
    const SplitQuaternion = makeSplitQuaternion(baseClass, 1, 2, 3);
    return SplitQuaternion;
  }
  if (p === 1 && q === 1 && r === 0) {
    let p1 = 1;
    let q_ = 2;
    if (baseClass.metric[0] < 0) {
      p1 = 2;
      q_ = 1;
    }

    const Coquaternion = makeSplitQuaternion(baseClass, p1, 3, q_);
    Object.defineProperty(Coquaternion, 'name', {value: 'Coquaternion'});
    return Coquaternion;
  }
  if (p === 0 && q === 2 && r === 0) {
    class Quaternion extends baseClass {
      get algebra() {
        return Quaternion;
      }

      sqrt(forceIter = false, numIter = 16): AlgebraElement {
        if (forceIter) {
          return super.sqrt(forceIter, numIter);
        }
        const result = this.imag();
        const imagNorm = result.vnorm();
        const [x, y] = complexSqrt(this.s, imagNorm);
        if (imagNorm < 1e-5) {
          result.rescale(0.5);
        } else {
          result.rescale(y / imagNorm);
        }
        result.s = x;
        return result;
      }

      exp(forceSeries = false, numTerms = 16): AlgebraElement {
        if (forceSeries) {
          return super.exp(forceSeries, numTerms);
        }
        const expS = Math.exp(this.s);
        const result = this.imag();
        const imagNorm = result.vnorm();
        result.rescale(expS * sinc(imagNorm));
        result.s = expS * Math.cos(imagNorm);
        return result;
      }

      log(): AlgebraElement {
        const imag = this.imag();
        const imagNorm = imag.vnorm();
        const result = imag.scale(Math.atan2(imagNorm, this.s) / imagNorm);
        result.s = Math.log(this.vnorm());
        return result;
      }

      inverse(): AlgebraElement {
        const conjugate = this.conjugate();
        return conjugate.scale(1 / this.mul(conjugate).s);
      }

      static zero() {
        return new Quaternion().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Quaternion(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Quaternion(baseClass.pseudoscalar(magnitude));
      }
      static basisBlade(...indices: number[]): AlgebraElement {
        return new Quaternion(baseClass.basisBlade(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Quaternion(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Quaternion(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Quaternion(baseClass.fromGanja(values));
      }
    }
    return Quaternion;
  }
  if (p === 1 && q === 0 && r === 1) {
    let r_ = 1;
    let pq = 2;
    if (baseClass.metric[0] !== 0) {
      pq = 1;
      r_ = 2;
    }
    const Euclidean1DPGA = makePGA1D(
      baseClass,
      splitComplexSqrt,
      Math.cosh,
      sinch,
      splitComplexLog,
      'Euclidean1DPGA',
      r_,
      pq
    );
    return Euclidean1DPGA;
  }
  if (p === 0 && q === 1 && r === 1) {
    let r_ = 1;
    let pq = 2;
    if (baseClass.metric[0] !== 0) {
      pq = 1;
      r_ = 2;
    }
    const ComplexPGA = makePGA1D(
      baseClass,
      complexSqrt,
      Math.cos,
      sinc,
      complexLog,
      'ComplexPGA',
      r_,
      pq
    );
    return ComplexPGA;
  }
  if (p === 0 && q === 0 && r === 2) {
    const DoublePGA = makePGA1D(
      baseClass,
      dualSqrt,
      y => 1, // eslint-disable-line @typescript-eslint/no-unused-vars
      y => 1, // eslint-disable-line @typescript-eslint/no-unused-vars
      dualLog,
      'DoublePGA'
    );
    return DoublePGA;
  }
  if (p === 2 && q === 0 && r === 1) {
    if (baseClass.metric[0] !== 0) {
      console.warn(
        'Custom metric order not supported for specialization in p=2 r=1'
      );
      return baseClass;
    }
    const SplitQuaternionPGA = makeSplitQuaternionPGA(baseClass, 2, 4, 6, 1);
    return SplitQuaternionPGA;
  }
  if (p === 1 && q === 1 && r === 1) {
    if (baseClass.metric[0] !== 0) {
      console.warn(
        'Custom metric order not supported for specialization in p=1 q=1 r=1'
      );
      return baseClass;
    }
    const CoquaternionPGA = makeSplitQuaternionPGA(baseClass, 2, 6, 4, -1);
    Object.defineProperty(CoquaternionPGA, 'name', {value: 'CoquaternionPGA'});
    return CoquaternionPGA;
  }
  if (p === 0 && q === 2 && r === 1) {
    if (baseClass.metric[0] !== 0) {
      console.warn(
        'Custom metric order not supported for specialization in q=2 r=1'
      );
      return baseClass;
    }
    class QuaternionPGA extends baseClass {
      get algebra() {
        return QuaternionPGA;
      }

      imagNorm() {
        return Math.sqrt(
          this[2] * this[2] + this[4] * this[4] + this[6] * this[6]
        );
      }

      sqrt(forceIter = false, numIter = 16): AlgebraElement {
        if (forceIter) {
          return super.sqrt(forceIter, numIter);
        }
        // Indices
        // 0 scalar
        // 1 e0
        // 2 i
        // 3 e0i
        // 4 j
        // 5 e0j
        // 6 ij = k
        // 7 e0ij = pseudoscalar

        const result = this.imag();
        const imagNorm = this.imagNorm();
        const [x, y] = complexSqrt(this.s, imagNorm);
        result.s = x;
        let im = 0.5;
        if (imagNorm > 1e-5) {
          im = y / imagNorm;
        }
        result[2] *= im;
        result[4] *= im;
        result[6] *= im;

        const ps =
          (-0.5 *
            (this[1] * result[6] +
              this[3] * result[4] -
              this[5] * result[2] -
              this.ps * x)) /
          (x * x + y * y);

        result[1] = (0.5 * this[1] + result[6] * ps) / x;
        result[3] = (0.5 * this[3] + result[4] * ps) / x;
        result[5] = (0.5 * this[5] - result[2] * ps) / x;
        result.ps = ps;

        return result;
      }

      inverse(): AlgebraElement {
        const reverse = this.rev();
        const involute = this.involute();
        const conjugate = this.conjugate();
        return reverse
          .mul(involute)
          .mul(conjugate)
          .scale(1 / this.mul(conjugate).mul(involute).mul(reverse).s);
      }

      static zero() {
        return new QuaternionPGA().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new QuaternionPGA(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new QuaternionPGA(baseClass.pseudoscalar(magnitude));
      }
      static basisBlade(...indices: number[]): AlgebraElement {
        return new QuaternionPGA(baseClass.basisBlade(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new QuaternionPGA(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new QuaternionPGA(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new QuaternionPGA(baseClass.fromGanja(values));
      }
    }
    return QuaternionPGA;
  }
  if (p === 3 && q === 0 && r === 0) {
    return make3D(baseClass, 1, 2, 4, 3, 5, 6, 'Elliptic2DPGA');
  }
  if (p === 1 && q === 2 && r === 0) {
    if (baseClass.metric[0] < 0) {
      console.warn(
        'Custom metric order not supported for specialization in p=1 q=2 r=0'
      );
      return baseClass;
    }
    return make3D(baseClass, 1, 5, 3, 4, 2, 6, 'Clifford3D');
  }
  if (p === 4 && q === 0 && r === 0) {
    class Elliptic3DPGA extends baseClass {
      get algebra() {
        return Elliptic3DPGA;
      }

      // Elliptic/Spherical PGA. e1*e1 = e2*e2 = e3*e3 = e4*e4 = 1
      // Normalize an even element X on the basis [1,e12,e13,e14,e23,e24,e34,e1234]
      rotorNormalize(): AlgebraElement {
        const X = this.rotor();

        const S =
          X[0] * X[0] +
          X[1] * X[1] +
          X[2] * X[2] +
          X[3] * X[3] +
          X[4] * X[4] +
          X[5] * X[5] +
          X[6] * X[6] +
          X[7] * X[7];
        const T = 2 * (X[0] * X[7] - X[1] * X[6] + X[2] * X[5] - X[3] * X[4]);
        const N = ((S * S - T * T) ** 0.5 + S) ** 0.5,
          N2 = N * N;
        const M = (2 ** 0.5 * N) / (N2 * N2 - T * T);
        const A = N2 * M,
          B = -T * M;
        return this.algebra.fromRotor([
          A * X[0] + B * X[7],
          A * X[1] - B * X[6],
          A * X[2] + B * X[5],
          A * X[3] - B * X[4],
          A * X[4] - B * X[3],
          A * X[5] + B * X[2],
          A * X[6] - B * X[1],
          A * X[7] + B * X[0],
        ]);
      }

      static zero() {
        return new Elliptic3DPGA().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Elliptic3DPGA(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Elliptic3DPGA(baseClass.pseudoscalar(magnitude));
      }
      static basisBlade(...indices: number[]): AlgebraElement {
        return new Elliptic3DPGA(baseClass.basisBlade(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Elliptic3DPGA(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Elliptic3DPGA(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Elliptic3DPGA(baseClass.fromGanja(values));
      }
    }
    return Elliptic3DPGA;
  }
  if (p === 3 && q === 1 && r === 0) {
    if (baseClass.metric[3] > 0) {
      console.warn(
        'Custom metric order not supported for specialization in p=3 q=1'
      );
      return baseClass;
    }
    class Hyperbolic3DPGA extends baseClass {
      get algebra() {
        return Hyperbolic3DPGA;
      }

      // STA/Hyperbolic PGA R3,1. e1*e1 = e2*e2 = e3*e3 = 1, e4*e4 = -1
      // Normalize an even element X on the basis [1,e12,e13,e14,e23,e24,e34,e1234]
      rotorNormalize(): AlgebraElement {
        const X = this.rotor();
        const S =
          X[0] * X[0] +
          X[1] * X[1] +
          X[2] * X[2] -
          X[3] * X[3] +
          X[4] * X[4] -
          X[5] * X[5] -
          X[6] * X[6] -
          X[7] * X[7];
        const T = 2 * (X[0] * X[7] - X[1] * X[6] + X[2] * X[5] - X[3] * X[4]);
        const N = ((S * S + T * T) ** 0.5 + S) ** 0.5,
          N2 = N * N;
        const M = (2 ** 0.5 * N) / (N2 * N2 + T * T);
        const A = N2 * M,
          B = -T * M;
        return this.algebra.fromRotor([
          A * X[0] - B * X[7],
          A * X[1] + B * X[6],
          A * X[2] - B * X[5],
          A * X[3] - B * X[4],
          A * X[4] + B * X[3],
          A * X[5] + B * X[2],
          A * X[6] - B * X[1],
          A * X[7] + B * X[0],
        ]);
      }

      static zero() {
        return new Hyperbolic3DPGA().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Hyperbolic3DPGA(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Hyperbolic3DPGA(baseClass.pseudoscalar(magnitude));
      }
      static basisBlade(...indices: number[]): AlgebraElement {
        return new Hyperbolic3DPGA(baseClass.basisBlade(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Hyperbolic3DPGA(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Hyperbolic3DPGA(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Hyperbolic3DPGA(baseClass.fromGanja(values));
      }
    }
    return Hyperbolic3DPGA;
  }
  if (p === 3 && q === 0 && r === 1) {
    if (baseClass.metric[0] !== 0) {
      console.warn(
        'Custom metric order not supported for specialization in p=3 r=1'
      );
      return baseClass;
    }
    class Euclidean3DPGA extends baseClass {
      get algebra() {
        return Euclidean3DPGA;
      }

      // 3D PGA. e1*e1 = e2*e2 = e3*e3 = 1, e0*e0 = 0
      // Normalize an even element X on the basis [1,e01,e02,e03,e12,e31,e23,e0123]
      rotorNormalize(): AlgebraElement {
        const X = this.rotor();
        const A =
          1 / (X[0] * X[0] + X[4] * X[4] + X[5] * X[5] + X[6] * X[6]) ** 0.5;
        const B =
          (X[7] * X[0] - (X[1] * X[6] + X[2] * X[5] + X[3] * X[4])) * A * A * A;
        return this.algebra.fromRotor([
          A * X[0],
          A * X[1] + B * X[6],
          A * X[2] + B * X[5],
          A * X[3] + B * X[4],
          A * X[4],
          A * X[5],
          A * X[6],
          A * X[7] - B * X[0],
        ]);
      }

      // Exponential of a bivector B (17 mul, 8 add, 2 div, 1 sincos, 1 sqrt)
      bivectorExp() {
        const B = this.vector(2);
        if (p === 3 && q === 0 && r === 1) {
          const l = B[3] * B[3] + B[4] * B[4] + B[5] * B[5];
          if (l === 0)
            return this.algebra.fromRotor([1, B[0], B[1], B[2], 0, 0, 0, 0]);
          const m = B[0] * B[5] + B[1] * B[4] + B[2] * B[3],
            a = Math.sqrt(l),
            c = Math.cos(a),
            s = Math.sin(a) / a,
            t = (m / l) * (c - s);
          return this.algebra.fromRotor([
            c,
            s * B[0] + t * B[5],
            s * B[1] + t * B[4],
            s * B[2] + t * B[3],
            s * B[3],
            s * B[4],
            s * B[5],
            m * s,
          ]);
        }
        return this.exp();
      }

      // Logarithm of a rotor R (14 mul, 5 add, 1 div, 1 acos, 1 sqrt)
      rotorLog() {
        const R = this.rotor();
        if (R[0] === 1)
          return this.algebra.fromVector([R[1], R[2], R[3], 0, 0, 0], 2);
        const a = 1 / (1 - R[0] * R[0]);
        const b = Math.acos(R[0]) * Math.sqrt(a);
        const c = a * R[7] * (1 - R[0] * b);
        return this.algebra.fromVector(
          [
            c * R[6] + b * R[1],
            c * R[5] + b * R[2],
            c * R[4] + b * R[3],
            b * R[4],
            b * R[5],
            b * R[6],
          ],
          2
        );
      }
      static zero() {
        return new Euclidean3DPGA().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Euclidean3DPGA(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Euclidean3DPGA(baseClass.pseudoscalar(magnitude));
      }
      static basisBlade(...indices: number[]): AlgebraElement {
        return new Euclidean3DPGA(baseClass.basisBlade(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Euclidean3DPGA(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Euclidean3DPGA(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Euclidean3DPGA(baseClass.fromGanja(values));
      }
    }
    return Euclidean3DPGA;
  }
  if (p === 4 && q === 1 && r === 0) {
    if (baseClass.metric[4] > 0) {
      console.warn(
        'Custom metric order not supported for specialization in p=4 q=1'
      );
      return baseClass;
    }
    class Conformal3DGA extends baseClass {
      get algebra() {
        return Conformal3DGA;
      }

      // CGA R4,1. e1*e1 = e2*e2 = e3*e3 = e4*4 = 1, e5*e5 = -1
      // Normalize an even element X = [1,e12,e13,e14,e15,e23,e24,e25,e34,e35,e45,e1234,e1235,e1245,e1345,e2345]
      rotorNormalize(): AlgebraElement {
        const X = this.rotor();

        const S =
          X[0] * X[0] -
          X[10] * X[10] +
          X[11] * X[11] -
          X[12] * X[12] -
          X[13] * X[13] -
          X[14] * X[14] -
          X[15] * X[15] +
          X[1] * X[1] +
          X[2] * X[2] +
          X[3] * X[3] -
          X[4] * X[4] +
          X[5] * X[5] +
          X[6] * X[6] -
          X[7] * X[7] +
          X[8] * X[8] -
          X[9] * X[9];
        const T1 =
          2 *
          (X[0] * X[11] -
            X[10] * X[12] +
            X[13] * X[9] -
            X[14] * X[7] +
            X[15] * X[4] -
            X[1] * X[8] +
            X[2] * X[6] -
            X[3] * X[5]);
        const T2 =
          2 *
          (X[0] * X[12] -
            X[10] * X[11] +
            X[13] * X[8] -
            X[14] * X[6] +
            X[15] * X[3] -
            X[1] * X[9] +
            X[2] * X[7] -
            X[4] * X[5]);
        const T3 =
          2 *
          (X[0] * X[13] -
            X[10] * X[1] +
            X[11] * X[9] -
            X[12] * X[8] +
            X[14] * X[5] -
            X[15] * X[2] +
            X[3] * X[7] -
            X[4] * X[6]);
        const T4 =
          2 *
          (X[0] * X[14] -
            X[10] * X[2] -
            X[11] * X[7] +
            X[12] * X[6] -
            X[13] * X[5] +
            X[15] * X[1] +
            X[3] * X[9] -
            X[4] * X[8]);
        const T5 =
          2 *
          (X[0] * X[15] -
            X[10] * X[5] +
            X[11] * X[4] -
            X[12] * X[3] +
            X[13] * X[2] -
            X[14] * X[1] +
            X[6] * X[9] -
            X[7] * X[8]);
        const TT = -T1 * T1 + T2 * T2 + T3 * T3 + T4 * T4 + T5 * T5;
        const N = ((S * S + TT) ** 0.5 + S) ** 0.5,
          N2 = N * N;
        const M = (2 ** 0.5 * N) / (N2 * N2 + TT);
        const A = N2 * M,
          [B1, B2, B3, B4, B5] = [-T1 * M, -T2 * M, -T3 * M, -T4 * M, -T5 * M];
        return this.algebra.fromRotor([
          A * X[0] +
            B1 * X[11] -
            B2 * X[12] -
            B3 * X[13] -
            B4 * X[14] -
            B5 * X[15],
          A * X[1] -
            B1 * X[8] +
            B2 * X[9] +
            B3 * X[10] -
            B4 * X[15] +
            B5 * X[14],
          A * X[2] +
            B1 * X[6] -
            B2 * X[7] +
            B3 * X[15] +
            B4 * X[10] -
            B5 * X[13],
          A * X[3] -
            B1 * X[5] -
            B2 * X[15] -
            B3 * X[7] -
            B4 * X[9] +
            B5 * X[12],
          A * X[4] -
            B1 * X[15] -
            B2 * X[5] -
            B3 * X[6] -
            B4 * X[8] +
            B5 * X[11],
          A * X[5] -
            B1 * X[3] +
            B2 * X[4] -
            B3 * X[14] +
            B4 * X[13] +
            B5 * X[10],
          A * X[6] +
            B1 * X[2] +
            B2 * X[14] +
            B3 * X[4] -
            B4 * X[12] -
            B5 * X[9],
          A * X[7] +
            B1 * X[14] +
            B2 * X[2] +
            B3 * X[3] -
            B4 * X[11] -
            B5 * X[8],
          A * X[8] -
            B1 * X[1] -
            B2 * X[13] +
            B3 * X[12] +
            B4 * X[4] +
            B5 * X[7],
          A * X[9] -
            B1 * X[13] -
            B2 * X[1] +
            B3 * X[11] +
            B4 * X[3] +
            B5 * X[6],
          A * X[10] +
            B1 * X[12] -
            B2 * X[11] -
            B3 * X[1] -
            B4 * X[2] -
            B5 * X[5],
          A * X[11] +
            B1 * X[0] +
            B2 * X[10] -
            B3 * X[9] +
            B4 * X[7] -
            B5 * X[4],
          A * X[12] +
            B1 * X[10] +
            B2 * X[0] -
            B3 * X[8] +
            B4 * X[6] -
            B5 * X[3],
          A * X[13] - B1 * X[9] + B2 * X[8] + B3 * X[0] - B4 * X[5] + B5 * X[2],
          A * X[14] + B1 * X[7] - B2 * X[6] + B3 * X[5] + B4 * X[0] - B5 * X[1],
          A * X[15] - B1 * X[4] + B2 * X[3] - B3 * X[2] + B4 * X[1] + B5 * X[0],
        ]);
      }

      static zero() {
        return new Conformal3DGA().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Conformal3DGA(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Conformal3DGA(baseClass.pseudoscalar(magnitude));
      }
      static basisBlade(...indices: number[]): AlgebraElement {
        return new Conformal3DGA(baseClass.basisBlade(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Conformal3DGA(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Conformal3DGA(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Conformal3DGA(baseClass.fromGanja(values));
      }
    }
    return Conformal3DGA;
  }

  return baseClass;
}
