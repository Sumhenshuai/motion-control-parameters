/* Writes tests/fixtures/cases.json: current-loop cases whose margins tools/gen_reference.py recomputes independently. */
const fs = require('fs'), path = require('path');
const core = require('../app/core.js');
const setups = [
  { name: 'firmware defaults (Renesas pole placement, Nd=1)', in: {} },
  { name: 'IMC 300 Hz, backward Euler, Nd=1', in: { curMethod: 'imc', fcc: '300', piImpl: 'inc', disc: 'be' } },
  { name: 'IMC 400 Hz, Tustin, Nd=1', in: { curMethod: 'imc', fcc: '400', piImpl: 'inc', disc: 'tu' } },
  { name: 'IMC 250 Hz, forward Euler, Nd=2', in: { curMethod: 'imc', fcc: '250', piImpl: 'inc', disc: 'fe', nd: '2' } },
  { name: 'z-domain cancellation 500 Hz, Nd=0', in: { curMethod: 'zd', fcc: '500', piImpl: 'inc', disc: 'be', nd: '0' } },
  { name: 'Renesas 150 Hz, double update, Nd=1', in: { fcc: '150', upd: 'double' } },
];
const cases = setups.map(s => {
  const R = core.computeAll(s.in);
  return ['d', 'q'].map(ax => ({
    name: `${s.name} · ${ax} axis`, inputs: s.in, axis: ax,
    Rs: R.P.Rs, L: R.cur[ax].L, Ts: R.P.Ts, Nd: R.P.Nd, Kp: R.cur[ax].Kp, Ki: R.cur[ax].Ki, b0: R.cur[ax].b0, b1: R.cur[ax].b1, disc: R.P.disc,
  }));
}).flat();
fs.writeFileSync(path.join(__dirname, '../tests/fixtures/cases.json'), JSON.stringify(cases, null, 1));
console.log(cases.length, 'cases');
