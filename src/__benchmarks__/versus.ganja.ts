import Benchmark = require("benchmark");
const GanjaAlgebra = require('ganja.js');
import Algebra, { AlgebraElement } from "..";


function randomElement(Ga: typeof AlgebraElement) {
  const value: number[] = [];
  while (value.length < Ga.size) {
    value.push(Math.random() * 4 - 2);
  }
  return new Ga(value);
}

const Cl3 = Algebra(3);
const a = randomElement(Cl3);
const b = randomElement(Cl3);

const GanjaCl3 = GanjaAlgebra(3);
const ganjaA = new GanjaCl3(a.ganja());
const ganjaB = new GanjaCl3(b.ganja());

const suite = new Benchmark.Suite;

suite.add('ts-geometric-algebra#mul', () => {
  a.mul(b);
})
.add('ganja.js#Mul', () => {
  ganjaA.Mul(ganjaB);
})
.on('cycle', (event: { target: any; }) => {
  console.log(String(event.target));
})
.on('complete', function() {
  console.log('Fastest is ' + suite.filter('fastest').map('name'));
})
.run({ 'async': true });
