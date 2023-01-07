export function copysign(x: number, y: number) {
  if (y > 0) {
    return Math.abs(x);
  }
  return -Math.abs(x);
}

export type Bivariate2D = (x: number, y: number) => [number, number];

export function complexSqrt(x: number, y: number): [number, number] {
  const r = Math.hypot(x, y);
  return [Math.sqrt((r + x) * 0.5), copysign(Math.sqrt((r - x) * 0.5), y)];
}

export function splitComplexSqrt(x: number, y: number): [number, number] {
  const r = Math.sqrt(x * x - y * y);
  return [Math.sqrt((x + r) * 0.5), copysign(Math.sqrt((x - r) * 0.5), y)];
}

export function dualSqrt(x: number, y: number): [number, number] {
  const s = Math.sqrt(x);
  return [s, (0.5 * y) / s];
}

export function complexExp(x: number, y: number): [number, number] {
  const expX = Math.exp(x);
  return [expX * Math.cos(y), expX * Math.sin(y)];
}

export function splitComplexExp(x: number, y: number): [number, number] {
  const expX = Math.exp(x);
  return [expX * Math.cosh(y), expX * Math.sinh(y)];
}

export function dualExp(x: number, y: number): [number, number] {
  const expX = Math.exp(x);
  return [expX, expX * y];
}

export function complexLog(x: number, y: number): [number, number] {
  const norm = Math.hypot(x, y);
  return [Math.log(norm), Math.atan2(y, x)];
}

export function splitComplexLog(x: number, y: number): [number, number] {
  const norm = Math.sqrt(x * x - y * y);
  return [Math.log(norm), Math.asinh(y / norm)];
}

export function dualLog(x: number, y: number): [number, number] {
  return [Math.log(x), y / x];
}

export function sinc(x: number) {
  if (Math.abs(x) < 1e-6) {
    return 1 - (x * x) / 6;
  }
  return Math.sin(x) / x;
}

export function sinch(x: number) {
  if (Math.abs(x) < 1e-6) {
    return 1 + (x * x) / 6;
  }
  return Math.sinh(x) / x;
}
