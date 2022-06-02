// Specialized implementations for particular metrics

import {AlgebraElement} from './element';
import {complexSqrt, splitComplexSqrt} from './utils';

// https://www.researchgate.net/publication/360528787_Normalization_Square_Roots_and_the_Exponential_and_Logarithmic_Maps_in_Geometric_Algebras_of_Less_than_6D
export function pqrMixin(
  p: number,
  q: number,
  r: number,
  baseClass: typeof AlgebraElement
): typeof AlgebraElement {
  if (p === 0 && q === 0 && r === 0) {
    class Scalar extends baseClass {
      cls() {
        return Scalar;
      }

      sqrt(forceBabylon = false, numIter = 16): AlgebraElement {
        if (forceBabylon) {
          return super.sqrt(forceBabylon, numIter);
        }
        return this.cls().scalar(Math.sqrt(this.s));
      }

      exp(forceTaylor = false, numTaylorTerms = 32): AlgebraElement {
        if (forceTaylor) {
          return super.exp(forceTaylor, numTaylorTerms);
        }
        return this.cls().scalar(Math.exp(this.s));
      }

      log(): AlgebraElement {
        return this.cls().scalar(Math.log(this.s));
      }

      inverse(): AlgebraElement {
        return this.cls().scalar(1 / this.s);
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
      static basisVector(...indices: number[]): AlgebraElement {
        return new Scalar(baseClass.basisVector(...indices));
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
    class SplitComplex extends baseClass {
      cls() {
        return SplitComplex;
      }

      sqrt(forceBabylon = false, numIter = 16): AlgebraElement {
        if (forceBabylon) {
          return super.sqrt(forceBabylon, numIter);
        }
        return new (this.cls())(splitComplexSqrt(this.s, this.ps));
      }

      exp(forceTaylor = false, numTaylorTerms = 32): AlgebraElement {
        if (forceTaylor) {
          return super.exp(forceTaylor, numTaylorTerms);
        }
        const expS = Math.exp(this.s);
        return new (this.cls())([
          expS * Math.cosh(this.ps),
          expS * Math.sinh(this.ps),
        ]);
      }

      log(): AlgebraElement {
        const norm = Math.sqrt(this.s ** 2 - this.ps ** 2);
        return new (this.cls())([Math.log(norm), Math.asinh(this.ps / norm)]);
      }

      inverse(): AlgebraElement {
        const involute = this.involute();
        return involute.scale(1 / this.mul(involute).s);
      }

      static zero() {
        return new SplitComplex().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new SplitComplex(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new SplitComplex(baseClass.pseudoscalar(magnitude));
      }
      static basisVector(...indices: number[]): AlgebraElement {
        return new SplitComplex(baseClass.basisVector(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new SplitComplex(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new SplitComplex(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new SplitComplex(baseClass.fromGanja(values));
      }
    }
    return SplitComplex;
  }
  if (p === 0 && q === 1 && r === 0) {
    class Complex extends baseClass {
      cls() {
        return Complex;
      }

      sqrt(forceBabylon = false, numIter = 16): AlgebraElement {
        if (forceBabylon) {
          return super.sqrt(forceBabylon, numIter);
        }
        return new (this.cls())(complexSqrt(this.s, this.ps));
      }

      exp(forceTaylor = false, numTaylorTerms = 32): AlgebraElement {
        if (forceTaylor) {
          return super.exp(forceTaylor, numTaylorTerms);
        }
        const expS = Math.exp(this.s);
        return new (this.cls())([
          expS * Math.cos(this.ps),
          expS * Math.sin(this.ps),
        ]);
      }

      log(): AlgebraElement {
        const norm = Math.hypot(this.s, this.ps);
        return new (this.cls())([Math.log(norm), Math.atan2(this.ps, this.s)]);
      }

      inverse(): AlgebraElement {
        const involute = this.involute();
        return involute.scale(1 / this.mul(involute).s);
      }

      static zero() {
        return new Complex().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Complex(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Complex(baseClass.pseudoscalar(magnitude));
      }
      static basisVector(...indices: number[]): AlgebraElement {
        return new Complex(baseClass.basisVector(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Complex(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Complex(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Complex(baseClass.fromGanja(values));
      }
    }
    return Complex;
  }
  if (p === 0 && q === 0 && r === 1) {
    class Dual extends baseClass {
      cls() {
        return Dual;
      }

      sqrt(forceBabylon = false, numIter = 16): AlgebraElement {
        if (forceBabylon) {
          return super.sqrt(forceBabylon, numIter);
        }
        const s = Math.sqrt(this.s);
        return new (this.cls())([s, (0.5 * this.ps) / s]);
      }

      exp(forceTaylor = false, numTaylorTerms = 32): AlgebraElement {
        if (forceTaylor) {
          return super.exp(forceTaylor, numTaylorTerms);
        }
        const expS = Math.exp(this.s);
        return new (this.cls())([expS, expS * this.ps]);
      }

      log(): AlgebraElement {
        return new (this.cls())([Math.log(this.s), this.ps / this.s]);
      }

      inverse(): AlgebraElement {
        const involute = this.involute();
        return involute.scale(1 / this.mul(involute).s);
      }

      static zero() {
        return new Dual().fill(0);
      }
      static scalar(magnitude = 1): AlgebraElement {
        return new Dual(baseClass.scalar(magnitude));
      }
      static pseudoscalar(magnitude = 1): AlgebraElement {
        return new Dual(baseClass.pseudoscalar(magnitude));
      }
      static basisVector(...indices: number[]): AlgebraElement {
        return new Dual(baseClass.basisVector(...indices));
      }
      static fromVector(
        values: Iterable<number>,
        grade?: number
      ): AlgebraElement {
        return new Dual(baseClass.fromVector(values, grade));
      }
      static fromRotor(values: Iterable<number>): AlgebraElement {
        return new Dual(baseClass.fromRotor(values));
      }
      static fromGanja(values: Iterable<number>) {
        return new Dual(baseClass.fromGanja(values));
      }
    }
    return Dual;
  }
  if (p === 0 && q === 2 && r === 0) {
    class Quaternion extends baseClass {
      cls() {
        return Quaternion;
      }

      sqrt(forceBabylon = false, numIter = 16): AlgebraElement {
        if (forceBabylon) {
          return super.sqrt(forceBabylon, numIter);
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

      exp(forceTaylor = false, numTaylorTerms = 32): AlgebraElement {
        if (forceTaylor) {
          return super.exp(forceTaylor, numTaylorTerms);
        }
        const expS = Math.exp(this.s);
        const result = this.imag();
        const imagNorm = result.vnorm();
        if (imagNorm < 1e-5) {
          result.rescale(expS);
        } else {
          result.rescale((expS * Math.sin(imagNorm)) / imagNorm);
        }
        result.s = expS * Math.cos(imagNorm);
        return result;
      }

      log(): AlgebraElement {
        const norm = this.vnorm();
        const imag = this.imag();
        const imagNorm = imag.vnorm();
        const result = imag.scale(Math.acos(this.s / norm) / imagNorm);
        result.s = Math.log(norm);
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
      static basisVector(...indices: number[]): AlgebraElement {
        return new Quaternion(baseClass.basisVector(...indices));
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
  if (p === 4 && q === 0 && r === 0) {
    class Elliptic3DPGA extends baseClass {
      cls() {
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
        return this.cls().fromRotor([
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
      static basisVector(...indices: number[]): AlgebraElement {
        return new Elliptic3DPGA(baseClass.basisVector(...indices));
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
    class Hyperbolic3DPGA extends baseClass {
      cls() {
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
        return this.cls().fromRotor([
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
      static basisVector(...indices: number[]): AlgebraElement {
        return new Hyperbolic3DPGA(baseClass.basisVector(...indices));
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
    class Euclidean3DPGA extends baseClass {
      cls() {
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
        return this.cls().fromRotor([
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
            return this.cls().fromRotor([1, B[0], B[1], B[2], 0, 0, 0, 0]);
          const m = B[0] * B[5] + B[1] * B[4] + B[2] * B[3],
            a = Math.sqrt(l),
            c = Math.cos(a),
            s = Math.sin(a) / a,
            t = (m / l) * (c - s);
          return this.cls().fromRotor([
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
          return this.cls().fromVector([R[1], R[2], R[3], 0, 0, 0], 2);
        const a = 1 / (1 - R[0] * R[0]);
        const b = Math.acos(R[0]) * Math.sqrt(a);
        const c = a * R[7] * (1 - R[0] * b);
        return this.cls().fromVector(
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
      static basisVector(...indices: number[]): AlgebraElement {
        return new Euclidean3DPGA(baseClass.basisVector(...indices));
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
    class Conformal3DGA extends baseClass {
      cls() {
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
        return this.cls().fromRotor([
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
      static basisVector(...indices: number[]): AlgebraElement {
        return new Conformal3DGA(baseClass.basisVector(...indices));
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
