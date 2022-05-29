import {describe, it, expect} from 'vitest';
import {decomposeQR, eigenValues, matrixCloseTo} from '../utils';

describe('QR decomposition', () => {
  it('decomposes a matrix to a known result', () => {
    const A = [
      [12, -51, 4],
      [6, 167, -68],
      [-4, 24, -41],
    ];

    const [Q, R] = decomposeQR(A);

    const expectedR = [
      [-14, -21, 14],
      [0, -175, 70],
      [0, 0, 35],
    ];
    const expectedQ = [
      [-6 / 7, 69 / 175, -58 / 175],
      [-3 / 7, -158 / 175, 6 / 175],
      [2 / 7, -6 / 35, -33 / 35],
    ];

    expect(matrixCloseTo(R, expectedR)).toBeTruthy();
    expect(matrixCloseTo(Q, expectedQ)).toBeTruthy();
  });
});

describe('Eigenvalue solver', () => {
  it('produces a known result', () => {
    const A = [
      [5, -10, -5],
      [2, 14, 2],
      [-4, -8, 6],
    ];

    const result = [10, 10, 5];

    eigenValues(A).forEach(eigenValue => {
      let i = 0;
      while (i < result.length) {
        if (Math.abs(eigenValue - result[i]) < 1e-5) {
          result.splice(i, 1);
        } else {
          ++i;
        }
      }
    });
    expect(result.length).toBe(0);
  });
});
