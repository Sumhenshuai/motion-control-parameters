/* Writes tests/fixtures/snapshot.json: key results for a few setups, so later changes that move a number are caught.
 * Regenerate deliberately (node tools/make_snapshot.js) only when a change in results is intended. */
const fs = require('fs'), path = require('path');
const core = require('../app/core.js');
const setups = {
  defaults: {},
  incremental: { piImpl: 'inc' },
  imc_bw: { curMethod: 'imc', fcc: '300', spdMethod: 'bw', fsc: '30', piImpl: 'inc', disc: 'tu' },
  spm_alt_hw: { curMethod: 'zd', fcc: '400', spdMethod: 'so', tf: 'amp', motorType: 'spm', Lq: '8.2', isnTopo: 'isnN', vdcTopo: 'vdcF', ocTopo: 'ocW', tmpTopo: 'tmpL' },
  fractional: { nd: '0.5', fcc: '150' },
};
const n = x => +(+x).toPrecision(8);
const snap = {};
for (const [k, s] of Object.entries(setups)) {
  const r = core.computeAll(s);
  snap[k] = { inputs: s, values: {
    curKp: n(r.cur.q.Kp), curKiTs: n(r.cur.q.Ki * r.P.Ts), curPM: n(r.cur.q.sw.pm), curBW: n(r.cur.q.sw.bw), curOS: n(r.cur.q.met.os),
    spdKp: n(r.spd.Kp), spdKi: n(r.spd.Ki), spdPM: n(r.spd.sw.pm), spdFc: n(r.spd.sw.fc), spdOS: n(r.spd.sim.met.os), spdDip: n(r.spd.sim.dip),
    posPM: n(r.pos.sw.pm), pllKp: n(r.obs.pll.Kp), emfKp: n(r.obs.emf.Kp), loL2: n(r.obs.lo.l2),
    hwGVA: n(r.hw.v.GVA), hwVoff: n(r.hw.v.Voff), hwKv: n(r.hw.v.Kv), hwTlsb: n(r.hw.v.tLsb), hwCnt: n(r.hw.v.cnt) } };
}
fs.writeFileSync(path.join(__dirname, '../tests/fixtures/snapshot.json'), JSON.stringify(snap, null, 1));
console.log('snapshot written:', Object.keys(snap).join(', '));
