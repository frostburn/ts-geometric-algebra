import {describe, it, expect} from 'vitest';

import Algebra, {AlgebraElement} from '../index';

function randomElement(Ga: typeof AlgebraElement) {
  const value: number[] = [];
  while (value.length < Ga.size) {
    value.push(Math.random() * 4 - 2);
  }
  return new Ga(value);
}

describe('TS-Geometric-Algebra unroll compatibility', () => {
  it('agrees with the unrolled methods using random elements', () => {
    for (let p = 0; p < 3; ++p) {
      for (let q = 0; q < 3; ++q) {
        for (let r = 0; r < 2; ++r) {
          const Ga = Algebra(p, q, r);
          const RolledGa = Algebra(p, q, r, Float32Array, false);

          for (let i = 0; i < 10; ++i) {
            const a = randomElement(Ga);
            const b = randomElement(Ga);

            const rA = new RolledGa(a);
            const rB = new RolledGa(b);

            expect(a.add(b).equals(rA.add(rB))).toBeTruthy();
            expect(a.sub(b).equals(rA.sub(rB))).toBeTruthy();
            expect(a.mul(b).closeTo(rA.mul(rB))).toBeTruthy();
            expect(a.wedge(b).closeTo(rA.wedge(rB))).toBeTruthy();
            expect(a.vee(b).closeTo(rA.vee(rB))).toBeTruthy();
            expect(a.dot(b).closeTo(rA.dot(rB))).toBeTruthy();
            expect(a.dotL(b).closeTo(rA.dotL(rB))).toBeTruthy();
          }
        }
      }
    }
  });
});
