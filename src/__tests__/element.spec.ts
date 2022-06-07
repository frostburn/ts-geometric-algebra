import {describe, it, expect} from 'vitest';
import {Algebra} from '..';
import {add, mul, linSolve} from '../element';

describe('Addition', () => {
  it('works on numbers', () => {
    const a = add(1, 2, 3);
    expect(a).toBe(6);
  });

  it('works on algebra elements', () => {
    const Grassmann = Algebra(0, 0, 2);
    const a = Grassmann.fromVector([1, 2]);
    const b = Grassmann.fromVector([2, 3]);
    const c = Grassmann.fromVector([4, 5]);
    const d = add(a, b, c);
    const expected = Grassmann.fromVector([7, 10]);
    expect(d.equals(expected)).toBeTruthy();
  });

  it('works on mixed elements', () => {
    const Grassmann = Algebra(0, 0, 2);
    const a = Grassmann.fromVector([1, 2]);
    const b = add(3, a, 4);
    const expected = new Grassmann([7, 1, 2, 0]);
    expect(b.equals(expected)).toBeTruthy();
  });
});

describe('Multiplication', () => {
  it('works on numbers', () => {
    const a = mul(1, 2, 3);
    expect(a).toBe(6);
  });

  it('works on algebra elements', () => {
    const Cl4 = Algebra(4);
    const a = Cl4.fromVector([1, -2, 3, 5]);
    const b = Cl4.fromVector([2, 3, -4, 6]);
    const c = Cl4.fromVector([4, 5, 6, -7]);
    const d = mul(a, b, c);
    const expected = Cl4.fromVector([59, 225, -137, -175]).add(
      Cl4.fromVector([88, -137, 246, 359], 3)
    );
    expect(d.equals(expected)).toBeTruthy();
  });

  it('works on mixed elements', () => {
    const Grassmann = Algebra(0, 0, 2);
    const a = Grassmann.fromVector([1, 2]);
    const b = mul(3, a, 4);
    const expected = Grassmann.fromVector([12, 24]);
    expect(b.equals(expected)).toBeTruthy();
  });
});

describe('Linear combination solver', () => {
  it('can solve a linear equation', () => {
    const Grassmann = Algebra(0, 0, 4);
    const a = Grassmann.fromVector([1, 2, 3, 4]);
    const b = Grassmann.fromVector([5, 4, 3, 2]);
    const c = Grassmann.fromVector([1, -1, 1, -1]);
    const d = Grassmann.fromVector([-2, 3, -1, 2]);

    const x = a.scale(2).add(b.scale(3)).add(c.scale(4)).add(d.scale(5));

    const coefficients = linSolve(x, [a, b, c, d]);
    expect(coefficients.length).toBe(4);
    expect(coefficients[0]).toBe(2);
    expect(coefficients[1]).toBe(3);
    expect(coefficients[2]).toBe(4);
    expect(coefficients[3]).toBe(5);
  });

  it('can solve a linear with a smaller basis than the embedding dimension', () => {
    const Grassmann = Algebra(0, 0, 5);
    const a = Grassmann.fromVector([-1, 2, -3, 4, 1]);
    const b = Grassmann.fromVector([5, 6, 7, 8, 9]);
    const c = Grassmann.fromVector([1, -1, 1, -1, 1]);

    const x = a.scale(-1).add(b.scale(3)).add(c.scale(-2));

    const coefficients = linSolve(x, [a, b, c]);
    expect(coefficients.length).toBe(3);
    expect(coefficients[0]).toBe(-1);
    expect(coefficients[1]).toBe(3);
    expect(coefficients[2]).toBe(-2);
  });

  it('returns NaNs when there is no solution', () => {
    const Grassmann = Algebra(0, 0, 4);
    const a = Grassmann.fromVector([-1, 2, -3, 4]);
    const b = Grassmann.fromVector([5, 6, 7, 8]);
    const c = Grassmann.fromVector([1, -1, 1, -1]);

    const x = Grassmann.fromVector([1, 2, 3, 4]);

    const coefficients = linSolve(x, [a, b, c]);
    expect(coefficients.length).toBe(3);
    expect(coefficients[0]).toBe(NaN);
    expect(coefficients[1]).toBe(NaN);
    expect(coefficients[2]).toBe(NaN);
  });
});
