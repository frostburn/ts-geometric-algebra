export function decomposeQR(matrix: number[][]) {
  const [m, n] = [matrix.length, matrix[0].length];
  const qr = matrix.map(row => [...row]);
  const Q = matrix.map(row => Array(row.length).fill(0));
  const R = matrix.map(row => Array(row.length).fill(0));
  const d = [];
  // helper matrix
  for (let k = 0; k < n; ++k) {
    let nrm = 0;
    for (let i = k; i < m; ++i) nrm += qr[i][k] * qr[i][k];
    if (nrm) {
      nrm = Math.sqrt(nrm);
      if (qr[k][k] < 0) nrm = -nrm;
      for (let i = k; i < m; ++i) qr[i][k] /= nrm;
      qr[k][k] = qr[k][k] + 1;
      for (let j = k + 1; j < n; ++j) {
        let s = 0;
        for (let i = k; i < m; ++i) s += qr[i][k] * qr[i][j];
        s = -s / qr[k][k];
        for (let i = k; i < m; ++i) qr[i][j] += s * qr[i][k];
      }
    }
    d.push(-nrm);
  }
  // extract Q
  for (let k = n - 1; k >= 0; --k) {
    for (let i = 0; i < m; ++i) Q[i][k] = 0;
    Q[k][k] = 1;
    if (qr[k][k]) {
      for (let j = k; j < n; ++j) {
        let s = 0;
        for (let i = k; i < m; ++i) s += qr[i][k] * Q[i][j];
        s = -s / qr[k][k];
        for (let i = k; i < m; ++i) Q[i][j] += s * qr[i][k];
      }
    }
  }
  // extract R
  for (let i = 0; i < n; ++i)
    for (let j = 0; j < n; ++j) R[i][j] = i < j ? qr[i][j] : i === j ? d[i] : 0;
  return [Q, R];
}

export function matrixCloseTo(A: number[][], B: number[][], tolerance = 1e-5) {
  if (A.length !== B.length) {
    return false;
  }
  for (let i = 0; i < A.length; ++i) {
    if (A[i].length !== B[i].length) {
      return false;
    }
    for (let j = 0; j < A[i].length; ++j) {
      if (Math.abs(A[i][j] - B[i][j]) > tolerance) {
        return false;
      }
    }
  }
  return true;
}

function mul(A: number[][], B: number[][]) {
  const res = A.map(r => Array(r.length).fill(0));
  for (let i = 0; i < A.length; ++i)
    for (let j = 0; j < B.length; ++j)
      for (let k = 0; k < B[0].length; ++k) res[i][k] += A[i][j] * B[j][k];
  return res;
}

export function eigenValues(matrix: number[][], numIter = 50) {
  for (let i = 0; i < numIter; ++i) {
    const [Q, R] = decomposeQR(matrix);
    matrix = mul(R, Q);
  }
  return matrix.map((x, i) => x[i]);
}

export function copysign(x: number, y: number) {
  if (y > 0) {
    return Math.abs(x);
  }
  return -Math.abs(x);
}

export function complexSqrt(x: number, y: number) {
  const r = Math.hypot(x, y);
  const sx = Math.sqrt((r + x) * 0.5);
  const sy = copysign(Math.sqrt((r - x) * 0.5), y);
  return [sx, sy];
}
