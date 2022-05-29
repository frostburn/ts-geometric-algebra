import {describe, it, expect} from 'vitest';

import Algebra, {AlgebraElement, ElementBaseType} from '../index';

function randomElement(Ga: typeof AlgebraElement) {
  const value: number[] = [];
  while (value.length < Ga.size) {
    value.push(Math.random() * 4 - 2);
  }
  return new Ga(value);
}

function randomVector(Ga: typeof AlgebraElement) {
  const value: number[] = [];
  while (value.length < Ga.dimensions) {
    value.push(Math.random() * 4 - 2);
  }
  return Ga.fromVector(value);
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

  it('fulfills the axioms of geometric algebra on random elements (positive metric)', () => {
    const dims = 4;
    const Cl4 = Algebra(dims);
    const scalar = Cl4.basisVector();
    const zero = Cl4.zero();

    function dot(a: ElementBaseType, b: ElementBaseType) {
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

      // Inner product
      const v = randomVector(Cl4);
      expect(
        v.mul(v).closeTo(scalar.scale(dot(v.vector(), v.vector())))
      ).toBeTruthy();
    }
  });

  it('fulfills the axioms of geometric algebra on random elements (mixed metric)', () => {
    const Ga = Algebra(1, 2, 1);
    const scalar = Ga.basisVector();
    const zero = Ga.zero();

    function dot(a: ElementBaseType, b: ElementBaseType) {
      return a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
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

      // Inner product
      const v = randomVector(Ga);
      expect(
        v.mul(v).closeTo(scalar.scale(dot(v.vector(), v.vector())))
      ).toBeTruthy();
    }
  });

  it('has an outer product that fulfills the axioms of geometric algebra on random elements', () => {
    const dims = 4;
    const Cl4 = Algebra(dims);
    const scalar = Cl4.basisVector();
    const zero = Cl4.zero();

    for (let i = 0; i < 10; ++i) {
      const elements: typeof scalar[] = [];
      for (let j = 0; j < 3; ++j) {
        elements.push(randomElement(Cl4));
      }
      const [a, b, c] = elements;

      // Closure
      expect(a.wedge(b).length).toBe(scalar.length);

      // Identity element
      expect(a.wedge(scalar).equals(a)).toBeTruthy();
      expect(scalar.wedge(a).equals(a)).toBeTruthy();

      // Associativity
      expect(a.wedge(b.wedge(c)).closeTo(a.wedge(b).wedge(c))).toBeTruthy();

      // Distributivity
      expect(
        a.wedge(b.add(c)).closeTo(a.wedge(b).add(a.wedge(c)))
      ).toBeTruthy();
      expect(
        b
          .add(c)
          .wedge(a)
          .closeTo(b.wedge(a).add(c.wedge(a)))
      ).toBeTruthy();

      const v = randomVector(Cl4);
      expect(v.wedge(v).closeTo(zero)).toBeTruthy();
    }
  });

  it('has a dual outer product that fulfills the axioms of geometric algebra on random elements', () => {
    const dims = 4;
    const Cl4 = Algebra(dims);
    const pseudoscalar = Cl4.pseudoscalar();
    const zero = Cl4.zero();

    for (let i = 0; i < 10; ++i) {
      const elements: typeof pseudoscalar[] = [];
      for (let j = 0; j < 3; ++j) {
        elements.push(randomElement(Cl4));
      }
      const [a, b, c] = elements;

      // Closure
      expect(a.vee(b).length).toBe(pseudoscalar.length);

      // Identity element
      expect(a.vee(pseudoscalar).equals(a)).toBeTruthy();
      expect(pseudoscalar.vee(a).equals(a)).toBeTruthy();

      // Associativity
      expect(a.vee(b.vee(c)).closeTo(a.vee(b).vee(c))).toBeTruthy();

      // Distributivity
      expect(a.vee(b.add(c)).closeTo(a.vee(b).add(a.vee(c)))).toBeTruthy();
      expect(
        b
          .add(c)
          .vee(a)
          .closeTo(b.vee(a).add(c.vee(a)))
      ).toBeTruthy();

      const v = randomVector(Cl4);
      expect(v.vee(v).closeTo(zero)).toBeTruthy();
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
        expect(element.div(element).closeTo(one)).toBeTruthy();
        expect(element.ldiv(element).closeTo(one)).toBeTruthy();
        expect(element.ldivs(element).closeTo(one)).toBeTruthy();
        break;
      }
    }
  });

  it('implements the hodge dual on basis vectors with a positive metric', () => {
    const Cl3 = Algebra(3);

    const scalar = Cl3.basisVector();

    const e0 = Cl3.basisVector(0);
    const e1 = Cl3.basisVector(1);
    const e2 = Cl3.basisVector(2);

    const e01 = Cl3.basisVector(0, 1);
    const e02 = Cl3.basisVector(0, 2);
    const e12 = Cl3.basisVector(1, 2);

    const pseudoscalar = Cl3.basisVector(0, 1, 2);

    expect(scalar.mul(scalar.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e0.mul(e0.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e1.mul(e1.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e2.mul(e2.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e01.mul(e01.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e02.mul(e02.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e12.mul(e12.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(
      pseudoscalar.mul(pseudoscalar.dual()).equals(pseudoscalar)
    ).toBeTruthy();

    expect(scalar.dual().equals(pseudoscalar)).toBeTruthy();
    expect(e0.dual().equals(e12)).toBeTruthy();
    expect(e1.dual().equals(e02.neg())).toBeTruthy();
    expect(e2.dual().equals(e01)).toBeTruthy();
    expect(e01.dual().equals(e2)).toBeTruthy();
    expect(e02.dual().equals(e1.neg())).toBeTruthy();
    expect(e12.dual().equals(e0)).toBeTruthy();
    expect(pseudoscalar.dual().equals(scalar)).toBeTruthy();
  });

  it('implements the hodge dual on basis vectors with a mixed metric', () => {
    const Cl3 = Algebra(1, 1, 1);

    const scalar = Cl3.basisVector();

    const e0 = Cl3.basisVector(0);
    const e1 = Cl3.basisVector(1);
    const e2 = Cl3.basisVector(2);

    const e01 = Cl3.basisVector(0, 1);
    const e02 = Cl3.basisVector(0, 2);
    const e12 = Cl3.basisVector(1, 2);

    const pseudoscalar = Cl3.basisVector(0, 1, 2);

    expect(scalar.mul(scalar.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e0.mul(e0.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e1.mul(e1.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e2.mul(e2.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e01.mul(e01.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e02.mul(e02.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(e12.mul(e12.dual()).equals(pseudoscalar)).toBeTruthy();
    expect(
      pseudoscalar.mul(pseudoscalar.dual()).equals(pseudoscalar)
    ).toBeTruthy();

    expect(scalar.dual().equals(pseudoscalar)).toBeTruthy();
    expect(e0.dual().equals(e12)).toBeTruthy();
    expect(e1.dual().equals(e02.neg())).toBeTruthy();
    expect(e2.dual().equals(e01)).toBeTruthy();
    expect(e01.dual().equals(e2)).toBeTruthy();
    expect(e02.dual().equals(e1.neg())).toBeTruthy();
    expect(e12.dual().equals(e0)).toBeTruthy();
    expect(pseudoscalar.dual().equals(scalar)).toBeTruthy();
  });

  it('satisfies operator definitions on random elements', () => {
    const dims = 6;
    const Cl4 = Algebra(dims);

    for (let i = 0; i < 10; ++i) {
      // ---Generic elements---
      const a = randomElement(Cl4);
      const b = randomElement(Cl4);

      // Unary
      expect(a.dual().undual().equals(a)).toBeTruthy();

      // Binary
      expect(a.mul(b).closeTo(b.rmul(a))).toBeTruthy();
      expect(a.wedge(b).closeTo(b.rwedge(a))).toBeTruthy();
      expect(a.vee(b).closeTo(b.dual().wedge(a.dual()).undual())).toBeTruthy();
      expect(a.vee(b).closeTo(b.rvee(a))).toBeTruthy();

      // ---Vector subspace---
      const u = randomVector(Cl4);
      const v = randomVector(Cl4);

      const inner = u.mul(v).add(v.mul(u)).scale(0.5);
      expect(u.dot(v).closeTo(inner)).toBeTruthy();
      expect(u.dotL(v).closeTo(inner)).toBeTruthy();
      expect(u.dotR(v).closeTo(inner)).toBeTruthy();
      expect(u.star(v).closeTo(inner)).toBeTruthy();
      expect(
        u.wedge(v).closeTo(u.mul(v).sub(v.mul(u)).scale(0.5))
      ).toBeTruthy();
    }
  });

  it('satisfies some identities on random elements', () => {
    const Ga = Algebra(3, 2);

    const pseudoscalar = Ga.pseudoscalar();

    for (let i = 0; i < 10; ++i) {
      const a = randomElement(Ga);
      const b = randomElement(Ga);
      const c = randomElement(Ga);

      expect(
        a.dotL(b).closeTo(a.wedge(b.div(pseudoscalar)).mul(pseudoscalar))
      ).toBeTruthy();
      expect(
        a
          .dotR(b)
          .closeTo(pseudoscalar.mul(pseudoscalar.inverse().mul(a).wedge(b)))
      ).toBeTruthy();
      expect(
        a
          .wedge(b)
          .star(c)
          .closeTo(a.star(b.dotL(c)))
      ).toBeTruthy();
      expect(c.star(b.wedge(a)).closeTo(c.dotR(b).star(a))).toBeTruthy();
      expect(a.dotL(b.dotL(c)).closeTo(a.wedge(b).dotL(c))).toBeTruthy();
      expect(
        a
          .dotL(b)
          .dotR(c)
          .closeTo(a.dotL(b.dotR(c)))
      ).toBeTruthy();
    }
  });

  it('can apply weights to the basis vectors', () => {
    const Cl3 = Algebra(3);
    const a = randomElement(Cl3);
    const b = a.applyWeights([0.5, 3, 5]);
    expect(a.s).toBe(b.s);
    const va = a.vector();
    const vb = b.vector();
    expect(vb[0]).toBeCloseTo(va[0] * 0.5);
    expect(vb[1]).toBeCloseTo(va[1] * 3);
    expect(vb[2]).toBeCloseTo(va[2] * 5);
    expect(b.ps).toBeCloseTo(a.ps * 0.5 * 3 * 5);
  });

  it('can read and write values in lexicographical order', () => {
    const Cl5 = Algebra(5);
    const scalar = Cl5.basisVector();
    const e = Array(5)
      .fill(null)
      .map((_, i) => Cl5.basisVector(i));

    const seven = Cl5.fromVector([7], 0);
    expect(seven.star(scalar).s).toBe(7);
    const bracketSeven = seven.vector(0);
    expect(bracketSeven.length).toBe(1);
    expect(bracketSeven[0]).toBe(7);

    const vector = Cl5.fromVector([1, 2, 3, 4, 5]);
    expect(vector.star(e[0]).s).toBe(1);
    expect(vector.star(e[1]).s).toBe(2);
    expect(vector.star(e[2]).s).toBe(3);
    expect(vector.star(e[3]).s).toBe(4);
    expect(vector.star(e[4]).s).toBe(5);
    const recovered = vector.vector();
    for (let i = 0; i < recovered.length; ++i) {
      expect(recovered[i]).toBe(i + 1);
    }

    const bivector = Cl5.fromVector(
      [-1, -2, -3, -4, -5, -6, -7, -8, -9, -10],
      2
    );
    expect(bivector.star(e[0].wedge(e[1])).s).toBe(1);
    expect(bivector.star(e[0].wedge(e[2])).s).toBe(2);
    expect(bivector.star(e[0].wedge(e[3])).s).toBe(3);
    expect(bivector.star(e[0].wedge(e[4])).s).toBe(4);
    expect(bivector.star(e[1].wedge(e[2])).s).toBe(5);
    expect(bivector.star(e[1].wedge(e[3])).s).toBe(6);
    expect(bivector.star(e[1].wedge(e[4])).s).toBe(7);
    expect(bivector.star(e[2].wedge(e[3])).s).toBe(8);
    expect(bivector.star(e[2].wedge(e[4])).s).toBe(9);
    expect(bivector.star(e[3].wedge(e[4])).s).toBe(10);
    const birecovered = bivector.vector(2);
    for (let i = 0; i < birecovered.length; ++i) {
      expect(birecovered[i]).toBe(-i - 1);
    }
  });

  it('can convert values in ganja order', () => {
    const Cl4 = Algebra(4);
    const element = randomElement(Cl4);
    const ganja = element.ganja();
    const recovered = Cl4.fromGanja(ganja);
    expect(element.equals(recovered)).toBeTruthy();
  });

  it('supports float64 precision', () => {
    const Cl3 = Algebra(1, 1, 1, Float64Array);
    const element = Cl3.zero();
    element.s = Number.MIN_VALUE;
    element.ps = Number.MAX_VALUE;
    expect(element.s).toBe(Number.MIN_VALUE);
    expect(element.ps).toBe(Number.MAX_VALUE);
  });

  it('implements the exponential function over the complex numbers', () => {
    const Complex = Algebra(0, 1);
    for (let i = 0; i < 10; ++i) {
      const z = randomElement(Complex);
      const expZ = z.exp(true);
      expect(expZ.s).toBeCloseTo(Math.exp(z.s) * Math.cos(z.ps));
      expect(expZ.ps).toBeCloseTo(Math.exp(z.s) * Math.sin(z.ps));
      const analytic = z.exp();
      expect(expZ.closeTo(analytic)).toBeTruthy();
    }
  });

  it('implements the exponential function over the quaternions', () => {
    const H = Algebra(0, 2);
    for (let i = 0; i < 10; ++i) {
      const z = randomElement(H);
      const expZ = z.exp(true);
      const imagZ = z.sub(z.grade(0));
      expect(expZ.s).toBeCloseTo(Math.exp(z.s) * Math.cos(imagZ.norm()));
      const imagExpZ = expZ.sub(expZ.grade(0));
      expect(
        imagExpZ.closeTo(
          imagZ.normalize(Math.exp(z.s) * Math.sin(imagZ.norm()))
        )
      ).toBeTruthy();
      const analytic = z.exp();
      expect(expZ.closeTo(analytic)).toBeTruthy();
    }
  });

  it('implements the exponential function over the hyperbolic numbers', () => {
    const Hyper = Algebra(1);
    for (let i = 0; i < 10; ++i) {
      const z = randomElement(Hyper);
      const expZ = z.exp(true);
      expect(expZ.s).toBeCloseTo(Math.exp(z.s) * Math.cosh(z.ps));
      expect(expZ.ps).toBeCloseTo(Math.exp(z.s) * Math.sinh(z.ps));
      const analytic = z.exp();
      expect(expZ.closeTo(analytic)).toBeTruthy();
    }
  });

  it('implements the exponential function over the dual numbers', () => {
    const Dual = Algebra(0, 0, 1);
    for (let i = 0; i < 10; ++i) {
      const z = randomElement(Dual);
      const expZ = z.exp(true);
      expect(expZ.s).toBeCloseTo(Math.exp(z.s));
      expect(expZ.ps).toBeCloseTo(Math.exp(z.s) * z.ps);
      const analytic = z.exp();
      expect(expZ.closeTo(analytic)).toBeTruthy();
    }
  });

  it('satisfies identities of the exponential function on random elements', () => {
    for (let p = 0; p < 3; ++p) {
      for (let q = 0; q < 3; ++q) {
        for (let r = 0; r < 2; ++r) {
          const Ga = Algebra(p, q, r);
          const a = randomElement(Ga);
          const b = randomElement(Ga);
          const c = a.scale(0.25).exp();
          const s = Ga.scalar(Math.random() * 4 - 2);

          expect(
            c.mul(b).div(c).exp().closeTo(c.mul(b.exp()).div(c), 0.1)
          ).toBeTruthy();
          // expect(a.mul(b).div(a).exp(true).closeTo(a.mul(b.exp(true)).div(a), 0.1)).toBeTruthy();

          expect(
            a.add(s).exp().closeTo(a.exp().mul(s.exp()), 0.1)
          ).toBeTruthy();

          const d = a.scale(0.25);
          const e = b.scale(0.25);
          expect(
            d
              .add(e)
              .exp()
              .closeTo(
                d
                  .scale(1 / 256)
                  .exp()
                  .mul(e.scale(1 / 256).exp())
                  .pow(256),
                0.1
              )
          ).toBeTruthy();
        }
      }
    }
  });

  it('can be raised to an integer power', () => {
    const Ga = Algebra(2, 1, 1);
    for (let i = 0; i < 10; ++i) {
      const a = randomElement(Ga);
      expect(a.pow(0).equals(Ga.scalar())).toBeTruthy();
      expect(a.pow(1).equals(a)).toBeTruthy();
      expect(a.pow(-1).equals(a.inverse())).toBeTruthy();
      expect(a.pow(2).equals(a.mul(a))).toBeTruthy();

      let power = a.mul(a).mul(a);
      expect(a.pow(3).closeTo(power, 0.01)).toBeTruthy();
      power = power.mul(a);
      expect(a.pow(4).closeTo(power, 0.05)).toBeTruthy();
      power = power.mul(a);
      expect(a.pow(5).closeTo(power, 0.1)).toBeTruthy();
      power = power.mul(a);
      expect(a.pow(6).closeTo(power, 0.1)).toBeTruthy();
      power = power.mul(a);
    }
  });

  it('can get and set individual incides', () => {
    const Ga = Algebra(3);
    const element = Ga.zero();
    expect(element.isNil()).toBeTruthy();
    element.setAt(1.25);
    expect(element.s).toBe(1.25);
    expect(element.isGrade(0)).toBeTruthy();
    element.s = 0;
    element.setAt(1, 0.5);
    expect(element.getAt(1)).toBe(0.5);
    expect(element.isGrade(1)).toBeTruthy();
    element.setAt(0, 2, 1.75);
    expect(element.getAt(0, 2)).toBe(1.75);
    expect(element.isNil()).toBeFalsy();
    expect(element.isGrade(0)).toBeFalsy();
    expect(element.isGrade(1)).toBeFalsy();
    expect(element.isGrade(2)).toBeFalsy();
    expect(element.isGrade(3)).toBeFalsy();
  });
});
