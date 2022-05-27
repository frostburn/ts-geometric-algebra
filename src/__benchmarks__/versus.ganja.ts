import Benchmark = require('benchmark');
const GanjaAlgebra = require('ganja.js');
import Algebra, {AlgebraElement} from '..';

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

const mulSuite = new Benchmark.Suite();
mulSuite
  .add('ts-geometric-algebra#mul', () => {
    a.mul(b);
  })
  .add('ganja.js#Mul', () => {
    ganjaA.Mul(ganjaB);
  })
  .on('cycle', (event: {target: any}) => {
    console.log(String(event.target));
  })
  .on('complete', () => {
    console.log('Fastest is ' + mulSuite.filter('fastest').map('name'));
  })
  .run({async: true});

const wedgeSuite = new Benchmark.Suite();
wedgeSuite
  .add('ts-geometric-algebra#wedge', () => {
    a.wedge(b);
  })
  .add('ganja.js#Wedge', () => {
    ganjaA.Wedge(ganjaB);
  })
  .on('cycle', (event: {target: any}) => {
    console.log(String(event.target));
  })
  .on('complete', () => {
    console.log('Fastest is ' + wedgeSuite.filter('fastest').map('name'));
  })
  .run({async: true});

const dotSuite = new Benchmark.Suite();
dotSuite
  .add('ts-geometric-algebra#dot', () => {
    a.dot(b);
  })
  .add('ganja.js#Dot', () => {
    ganjaA.Dot(ganjaB);
  })
  .on('cycle', (event: {target: any}) => {
    console.log(String(event.target));
  })
  .on('complete', () => {
    console.log('Fastest is ' + dotSuite.filter('fastest').map('name'));
  })
  .run({async: true});

const veeSuite = new Benchmark.Suite();
veeSuite
  .add('ts-geometric-algebra#vee', () => {
    a.vee(b);
  })
  .add('ganja.js#Vee', () => {
    ganjaA.Vee(ganjaB);
  })
  .on('cycle', (event: {target: any}) => {
    console.log(String(event.target));
  })
  .on('complete', () => {
    console.log('Fastest is ' + veeSuite.filter('fastest').map('name'));
  })
  .run({async: true});

const inverseSuite = new Benchmark.Suite();
inverseSuite
  .add('ts-geometric-algebra#inverse', () => {
    a.inverse();
  })
  .add('ganja.js#Inverse', () => {
    ganjaA.Inverse;
  })
  .on('cycle', (event: {target: any}) => {
    console.log(String(event.target));
  })
  .on('complete', () => {
    console.log('Fastest is ' + inverseSuite.filter('fastest').map('name'));
  })
  .run({async: true});
