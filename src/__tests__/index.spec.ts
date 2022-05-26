import {describe, it, expect} from 'vitest';

import Algebra, {AlgebraElement} from '../index';

function randomElement(Ga: typeof AlgebraElement) {
  const value: number[] = [];
  while (value.length < Ga.size) {
    value.push(Math.random() * 4 - 2);
  }
  return new Ga(value);
}

describe('Geometric Algebra', () => {
  it('supports a positive metric', () => {
    const Cl3 = Algebra(3);

    const scalar = Cl3.basisVector();

    const e0 = Cl3.basisVector(0);
    const e1 = Cl3.basisVector(1);
    const e2 = Cl3.basisVector(2);

    const e01 = Cl3.basisVector(0, 1);
    const e02 = Cl3.basisVector(0, 2);
    const e12 = Cl3.basisVector(1, 2);

    const pseudoscalar = Cl3.basisVector(0, 1, 2);

    expect(e0.mul(scalar).equals(e0)).toBeTruthy();
    expect(e0.mul(e0).equals(scalar)).toBeTruthy();
    expect(e0.mul(e1).equals(e01)).toBeTruthy();
    expect(e0.mul(e2).equals(e02)).toBeTruthy();
    expect(e0.mul(e01).equals(e1)).toBeTruthy();
    expect(e0.mul(e02).equals(e2)).toBeTruthy();
    expect(e0.mul(e12).equals(pseudoscalar)).toBeTruthy();
    expect(e0.mul(pseudoscalar).equals(e12)).toBeTruthy();

    expect(e1.mul(e0).equals(e01.neg())).toBeTruthy();
    expect(e1.mul(e1).equals(scalar)).toBeTruthy();
    expect(e1.mul(e2).equals(e12)).toBeTruthy();
    expect(e1.mul(e01).equals(e0.neg())).toBeTruthy();
    expect(e1.mul(e02).equals(pseudoscalar.neg())).toBeTruthy();
    expect(e1.mul(e12).equals(e2)).toBeTruthy();
    expect(e1.mul(pseudoscalar).equals(e02.neg())).toBeTruthy();

    expect(e2.mul(e0).equals(e02.neg())).toBeTruthy();
    expect(e2.mul(e1).equals(e12.neg())).toBeTruthy();
    expect(e2.mul(e2).equals(scalar)).toBeTruthy();

    expect(pseudoscalar.mul(pseudoscalar).equals(scalar.neg())).toBeTruthy();
  });

  it('fulfils the axioms of geometric algebra on random elements (positive metric)', () => {
    const dims = 4;
    const Cl4 = Algebra(dims);
    const scalar = Cl4.basisVector();
    const zero = Cl4.zero();

    function dot(a: number[], b: number[]) {
      let result = 0;
      for (let i = 0; i < dims; ++i) {
        result += a[i] * b[i];
      }
      return result;
    }

    // At least two distinct elements
    expect(!scalar.equals(zero)).toBeTruthy();

    for (let i = 0; i < 10; ++i) {
      const elements: typeof scalar[] = [];
      for (let j = 0; j < 3; ++j) {
        elements.push(randomElement(Cl4));
      }
      const [a, b, c] = elements;

      // Closure
      expect(a.mul(b).length).toBe(scalar.length);

      // Identity element
      expect(a.mul(scalar).equals(a)).toBeTruthy();
      expect(scalar.mul(a).equals(a)).toBeTruthy();

      // Associativity
      expect(a.mul(b.mul(c)).closeTo(a.mul(b).mul(c))).toBeTruthy();

      // Distributivity
      expect(a.mul(b.add(c)).closeTo(a.mul(b).add(a.mul(c)))).toBeTruthy();
      expect(
        b
          .add(c)
          .mul(a)
          .closeTo(b.mul(a).add(c.mul(a)))
      ).toBeTruthy();

      const value: number[] = [];
      while (value.length < dims) {
        value.push(Math.random() * 4 - 2);
      }
      const v = Cl4.fromVector(value);
      expect(v.mul(v).closeTo(scalar.scale(dot(value, value)))).toBeTruthy();
    }
  });

  it('fulfils the axioms of geometric algebra on random elements (mixed metric)', () => {
    const dims = 4;
    const Ga = Algebra(1, 2, 1);
    const scalar = Ga.basisVector();
    const zero = Ga.zero();

    function dot(a: number[], b: number[]) {
      return a[0] * b[0] - a[1] * b[1] - a[2] * b[2];
    }

    // At least two distinct elements
    expect(!scalar.equals(zero)).toBeTruthy();

    for (let i = 0; i < 10; ++i) {
      const elements: typeof scalar[] = [];
      for (let j = 0; j < 3; ++j) {
        elements.push(randomElement(Ga));
      }
      const [a, b, c] = elements;

      // Closure
      expect(a.mul(b).length).toBe(scalar.length);

      // Identity element
      expect(a.mul(scalar).equals(a)).toBeTruthy();
      expect(scalar.mul(a).equals(a)).toBeTruthy();

      // Associativity
      expect(a.mul(b.mul(c)).closeTo(a.mul(b).mul(c))).toBeTruthy();

      // Distributivity
      expect(a.mul(b.add(c)).closeTo(a.mul(b).add(a.mul(c)))).toBeTruthy();
      expect(
        b
          .add(c)
          .mul(a)
          .closeTo(b.mul(a).add(c.mul(a)))
      ).toBeTruthy();

      const value: number[] = [];
      while (value.length < dims) {
        value.push(Math.random() * 4 - 2);
      }
      const v = Ga.fromVector(value);
      expect(v.mul(v).closeTo(scalar.scale(dot(value, value)))).toBeTruthy();
    }
  });

  it('implements inverse in 0 dimensions', () => {
    const Cl0 = Algebra(0);
    const one = Cl0.scalar();
    const element = Cl0.scalar(Math.random() + 4);
    expect(element.mul(element.inverse()).closeTo(one)).toBeTruthy();
  });

  it('implements inverse in 1 to 6 dimensions', () => {
    for (let n = 1; n <= 6; ++n) {
      const Cl = Algebra(n);
      const zero = Cl.zero();
      const one = Cl.scalar();
      while (true) {
        const element = randomElement(Cl);
        const inverse = element.inverse();
        if (
          inverse.hasNaN() ||
          inverse.hasInfinity() ||
          inverse.closeTo(zero)
        ) {
          continue;
        }
        expect(element.mul(inverse).closeTo(one)).toBeTruthy();
        break;
      }
    }
  });

  it('implements the hodge dual on basis vectors', () => {
    const Cl3 = Algebra(3);

    const scalar = Cl3.basisVector();

    const e0 = Cl3.basisVector(0);
    const e1 = Cl3.basisVector(1);
    const e2 = Cl3.basisVector(2);

    const e01 = Cl3.basisVector(0, 1);
    const e02 = Cl3.basisVector(0, 2);
    const e12 = Cl3.basisVector(1, 2);

    const pseudoscalar = Cl3.basisVector(0, 1, 2);

    expect(scalar.dual().equals(pseudoscalar)).toBeTruthy();
    expect(e0.dual().equals(e12)).toBeTruthy();
    expect(e1.dual().equals(e02.neg())).toBeTruthy();
    expect(e2.dual().equals(e01)).toBeTruthy();
    expect(e01.dual().equals(e2.neg())).toBeTruthy();
    expect(e02.dual().equals(e1)).toBeTruthy();
    expect(e12.dual().equals(e0.neg())).toBeTruthy();
    expect(pseudoscalar.dual().equals(scalar.neg())).toBeTruthy();
  });

  it('satisfies operator definitions on random elements', () => {
    const dims = 4;
    const Cl4 = Algebra(dims);
    const pseudoscalar = Cl4.pseudoscalar();

    for (let i = 0; i < 10; ++i) {
      const a = randomElement(Cl4);
      const b = randomElement(Cl4);

      // Unary
      expect(a.dual().equals(a.mul(pseudoscalar))).toBeTruthy();
      expect(a.dual().undual().equals(a)).toBeTruthy();

      // Binary
      expect(a.dot(b).closeTo(a.mul(b).add(b.mul(a)).scale(0.5))).toBeTruthy();
      expect(
        a.wedge(b).closeTo(a.mul(b).sub(b.mul(a)).scale(0.5))
      ).toBeTruthy();
      expect(a.vee(b).closeTo(a.dual().wedge(b.dual()).dual())).toBeTruthy();
    }
  });
});
