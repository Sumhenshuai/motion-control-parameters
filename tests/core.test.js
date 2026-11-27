/* Tests for app/core.js. Run: node --test tests/   (Node 18+, no dependencies) */
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('fs');
const path = require('path');
const core = require('../app/core.js');

const fx = f => path.join(__dirname, 'fixtures', f);
const close = (a, b, rel, msg) => {
  const tol = Math.max(rel * Math.abs(b), 1e-12);
  assert.ok(Math.abs(a - b) <= tol, `${msg}: got ${a}, expected ${b} (±${tol})`);
};
const TAU = 2 * Math.PI;

test('PI discretisation coefficients', () => {
  const Kp = 2, Ki = 500, T = 1e-4;
  assert.deepEqual(core.piCoef(Kp, Ki, T, 'be'), { b0: Kp + Ki * T, b1: -Kp });
  assert.deepEqual(core.piCoef(Kp, Ki, T, 'fe'), { b0: Kp, b1: -Kp + Ki * T });
  const tu = core.piCoef(Kp, Ki, T, 'tu');
  close(tu.b0, Kp + Ki * T / 2, 1e-15, 'tustin b0'); close(tu.b1, -Kp + Ki * T / 2, 1e-15, 'tustin b1');
});

test('current-loop margins and poles match python-control (independent reference)', () => {
  const ref = JSON.parse(fs.readFileSync(fx('reference.json'), 'utf8'));
  assert.ok(ref.length >= 10);
  for (const r of ref) {
    const R = core.computeAll(r.inputs), a = R.cur[r.axis];
    close(a.Kp, r.Kp, 1e-12, `${r.name} Kp`);
    assert.ok(Math.abs(a.sw.pm - r.pm) < 0.02, `${r.name} PM ${a.sw.pm} vs ${r.pm}`);
    close(a.sw.fc, r.fc_hz, 1e-3, `${r.name} crossover`);
    if (r.gm_db === null) assert.ok(!isFinite(a.sw.gm) || a.sw.gm > 40, `${r.name} GM should be infinite, got ${a.sw.gm}`);
    else assert.ok(Math.abs(a.sw.gm - r.gm_db) < 0.05, `${r.name} GM ${a.sw.gm} vs ${r.gm_db}`);
    const maxPole = Math.max(...a.poles.map(p => Math.hypot(p.re, p.im)));
    close(maxPole, r.max_pole, 1e-9, `${r.name} max |pole|`);
  }
});

/* The time-domain simulation integrates the continuous RL plant in sub-steps and applies each sampled voltage
 * Nd·Ts after sampling. The z-domain model uses the modified z-transform. Both are derived separately, so their
 * sample-instant responses agreeing checks the fractional-delay model. */
for (const nd of ['0', '0.25', '0.5', '1', '1.5', '2']) {
  for (const impl of ['inc', 'rn']) {
    test(`fractional delay: z-model equals time-domain simulation (Nd=${nd}, ${impl})`, () => {
      const R = core.computeAll({ nd, piImpl: impl, curMethod: 'imc', fcc: '200', disc: 'be', curStep: '0.05' });
      const P = R.P, ax = R.cur.d, pl = core.plantZ(P, ax.L);
      close(pl.g0 + pl.g1, (1 - pl.a) / P.Rs, 1e-12, 'g0+g1 equals ZOH gain');
      // closed loop T = N / (D + N) with N = (b0 z + b1)·num, D = (z − 1)·den
      const mul = (a, b) => { const r = new Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => { r[i + j] += x * y; })); return r; };
      const N = mul([ax.b0, ax.b1], pl.num), D0 = mul([1, -1], pl.den);
      const n = D0.length - 1, Np = Array(n + 1 - N.length).fill(0).concat(N);
      const D = D0.map((d, i) => d + Np[i]);
      const amp = ax.sim.amp, K = Math.min(60, Math.floor(ax.sim.y.length / 5)), y = [], u = k => (k >= 0 ? amp : 0);
      for (let k = 0; k < K; k++) {
        let acc = 0;
        for (let j = 0; j <= n; j++) acc += Np[j] * u(k - j);
        for (let j = 1; j <= n; j++) acc -= D[j] * (k - j >= 0 ? y[k - j] : 0);
        y.push(acc / D[0]);
      }
      // simulation stores 5 points per sample period; index 5k+4 is the sample instant (k+1)·Ts
      assert.ok(K >= 30, 'simulation long enough');
      for (let k = 0; k < K - 1; k++) {
        const ySim = ax.sim.y[5 * k + 4], yZ = y[k + 1];
        assert.ok(Math.abs(ySim - yZ) < 1e-9 * amp + 1e-12, `k=${k}: sim ${ySim} vs z ${yZ}`);
      }
    });
  }
}

test('Renesas firmware gain formulas (TG55L, CURRENT_CFG_OMEGA 20 Hz, SPEED_CFG_OMEGA 2 Hz)', () => {
  const R = core.computeAll({});
  const w = TAU * 20, Ts = 1 / 5000;
  close(R.cur.d.Kp, 2 * 1 * w * 0.0082 - 1.34, 1e-12, 'kp_d');
  close(R.cur.q.Kp, 2 * 1 * w * 0.01123 - 1.34, 1e-12, 'kp_q');
  close(R.cur.q.Ki * Ts, w * w * 0.01123 * Ts, 1e-12, 'ki_q (firmware stores Ki·Ts)');
  // speed: firmware uses Kt = p·ψa with equal-power currents; the tool converts with kT = √1.5
  const ws = TAU * 2, Tw = 5 * Ts, ktFw = 3 * 0.086;
  close(R.spd.Kp * R.P.kT, 2 * ws * 3e-4 / ktFw, 1e-9, 'speed kp');
  close(R.spd.Ki * R.P.kT * Tw, ws * ws * 3e-4 / ktFw * Tw, 1e-9, 'speed ki');
  // Renesas cfg export reproduces the firmware macros
  const cfg = core.renesasCfg(R);
  assert.match(cfg, /CURRENT_CFG_OMEGA\s+\(20\.0f\)/);
  assert.match(cfg, /SPEED_CFG_OMEGA\s+\(2\.0f\)/);
  assert.match(cfg, /MOTOR_CFG_MAGNETIC_FLUX\s+\(0\.086f\)/);
});

test('PI implementations: Renesas form shares one limit between integral and output', () => {
  const pi = core.makePI('rn', 1, 0.5, 0, 0, 2);
  assert.equal(pi.step(10), 2);          // output saturates
  for (let i = 0; i < 10; i++) pi.step(10);
  assert.equal(pi.step(-1), 0.5);        // integral held at 2, then 2 − 0.5 = 1.5 → output −1 + 1.5
  const inc = core.makePI('inc', 0, 0, 1.5, -1, 3);
  assert.equal(inc.step(1), 1.5); assert.equal(inc.step(1), 2); assert.equal(inc.step(10), 3);
});

test('hardware conversions for the ARISCO V1.2 board', () => {
  const v = core.computeAll({}).hw.v;
  close(v.gainFw, 10000 / 2210, 1e-12, 'differential gain Rf/(Rg+Rser)');
  close(v.Voff, 2.5, 1e-12, 'bias');
  close(v.offCode, 2047.5, 1e-12, 'zero-current code');
  close(v.aLsb, (5 / 4095) / (10000 / 2210 * 0.02), 1e-12, 'A per LSB');
  close(v.Kv, 100.6, 1e-12, 'bus divider');
  close(v.cnt, 8000, 1e-12, 'PWM half-period count');
  const oc = core.computeAll({}).hw.out.oc.find(r => r[0] === '比较器动作电流');
  close(oc[2], 5 * 4.99 / 51.99 / 0.015, 1e-12, 'comparator trip current');
});

test('tolerance analysis: shunt-only tolerance gives exactly that gain error', () => {
  const zero = { tolR: '0', tolSh: '0', tolV: '0', tolSc: '0', tolNtc: '0', tolB: '0', opVos: '0' };
  let v = core.computeAll(zero).hw.v;
  close(v.gDev, 0, 1e-12, 'no tolerance → no gain error');
  v = core.computeAll({ ...zero, tolSh: '1' }).hw.v;
  close(v.gDev, 1, 1e-9, 'gain ∝ shunt');
  const oc = core.computeAll({ ...zero, tolSh: '1' }).hw.out.oc;
  const nom = 5 * 4.99 / 51.99 / 0.015;
  close(oc.find(r => r[0] === '比较器动作电流 最小 (容差)')[2], nom / 1.01, 1e-12, 'min trip with +1% shunt');
  close(oc.find(r => r[0] === '比较器动作电流 最大 (容差)')[2], nom / 0.99, 1e-12, 'max trip with −1% shunt');
  // op-amp offset alone: zero shift = Vos × noise gain
  const isn = core.computeAll({ ...zero, opVos: '2' }).hw.out.isn;
  close(isn.find(r => r[0] === '零点最坏偏移')[2], 2 * (1 + 10000 / 2210), 1e-9, 'Vos × (1 + Rf/Rn) in mV');
});

test('PLL speed estimate: unity at low frequency, matches continuous model well below fs', () => {
  const R = core.computeAll({});
  const g = R.spd.pll, P = R.P;
  const lo = core.pllResp(P, g, TAU * 0.01);
  close(Math.hypot(lo.re, lo.im), 1, 1e-4, '|ω̂/ω| at 0.01 Hz');
  for (const f of [5, 20, 40]) {
    const s = { re: 0, im: TAU * f }, num = { re: g.Ki, im: g.Kp * s.im }, den = { re: g.Ki - s.im * s.im, im: g.Kp * s.im };
    const hc = Math.hypot(num.re, num.im) / Math.hypot(den.re, den.im);
    const hd = core.pllResp(P, g, TAU * f);
    close(Math.hypot(hd.re, hd.im), hc, 0.02, `|ω̂/ω| at ${f} Hz vs continuous`);
  }
});

test('auto-tune lands on the phase-margin target and is the limit', () => {
  for (const [loop, pm] of [['cur', 60], ['spd', 50]]) {
    const r = core.autoTune({}, loop, pm, 6);
    assert.ok(r && r.value > 0, `${loop} tuned`);
    assert.ok(Math.abs(r.pm - pm) < 0.05, `${loop} PM ${r.pm}`);
    assert.ok(r.gm >= 6);
    const over = core.computeAll({ [r.key]: String(r.value * 1.05) });
    const pmOver = loop === 'cur' ? Math.min(over.cur.d.sw.pm, over.cur.q.sw.pm) : over.spd.sw.pm;
    assert.ok(pmOver < pm, `${loop}: 5% more bandwidth must break the target (${pmOver})`);
  }
});

test('worst-case auto-tune meets the target at every plant corner', () => {
  for (const [loop, pm] of [['cur', 60], ['spd', 50]]) {
    const r = core.autoTune({}, loop, pm, 6, 'worst');
    assert.ok(r && r.value > 0 && Math.abs(r.worstPm - pm) < 0.1, `${loop} worst PM ${r && r.worstPm}`);
    const nom = core.autoTune({}, loop, pm, 6, 'nom');
    assert.ok(r.value <= nom.value, 'robust tuning cannot exceed the nominal limit');
    const rb = core.robustness({ [r.key]: String(r.value) }, loop);
    assert.ok(rb.worst.pm >= pm - 0.5, `robustness table agrees (${rb.worst.pm})`);
  }
});

test('robustness: nominal corner equals the designed loop; zero ranges collapse to nominal', () => {
  const R = core.computeAll({});
  const rc = core.robustness({}, 'cur');
  assert.ok(Math.abs(rc.nominal.pm - Math.min(R.cur.d.sw.pm, R.cur.q.sw.pm)) < 0.05);
  assert.ok(rc.worst.pm <= rc.nominal.pm);
  const rs0 = core.robustness({ rbL: '0', rbKt: '0', rbJlo: '1', rbJhi: '1', rbRlo: '1', rbRhi: '1' }, 'spd');
  assert.equal(rs0.rows.length, 1);
  assert.ok(Math.abs(rs0.worst.pm - R.spd.sw.pm) < 0.2, `${rs0.worst.pm} vs ${R.spd.sw.pm}`);
});

test('NTC divider (B model)', () => {
  const R = core.computeAll({ tmpTopo: 'tmpL', tot: '100' });
  const row = R.hw.out.tmp.find(r => r[0] === '过温阈值电压');
  const Rn = 10000 * Math.exp(3950 * (1 / 373.15 - 1 / 298.15));
  close(row[2], 5 * Rn / (10000 + Rn), 1e-12, 'V at 100 °C');
});

test('C-style expression evaluator', () => {
  assert.equal(core.cArith('1 + 2 * 3'), 7);
  assert.equal(core.cArith('(1 + 2) * 3 / 2'), 4.5);
  assert.equal(core.cArith('-(2) * -3'), 6);
  assert.equal(core.cArith('1 && (0 || 1)'), 1);
  assert.equal(core.cArith('!1 || 2 == 2'), 1);
  assert.equal(core.cArith('3 >= 4'), 0);
  assert.ok(Number.isNaN(core.cArith('1 +')));
  assert.ok(Number.isNaN(core.cArith('alert(1)')));
});

test('cfg header parser: comments, conditionals, suffixes, hex and expressions', () => {
  const { vals, unresolved } = core.parseCfg(fs.readFileSync(fx('sample_cfg.h'), 'utf8'));
  assert.equal(vals.MOTOR_CFG_POLE_PAIRS, 4);                 // #ifdef branch taken, #else ignored
  assert.equal(vals.MOTOR_MCU_CFG_CARRIER_FREQ, 8);           // #elif branch (1-shunt define is commented out)
  assert.equal(vals.MOTOR_MCU_CFG_ADC_OFFSET, 2047);          // hex
  assert.equal(vals.MOTOR_COMMON_SPEED_LOOP_DIVIDER, 4);      // U suffix
  assert.equal(vals.CURRENT_CFG_PI_INTEGRAL_LIMIT_VQ, 150);   // expression with another macro
  close(vals.MOTOR_COMMON_PI_INTEGRAL_LIMIT_IQ, 6 * 1.7320508, 1e-9, 'expression with constant');
  close(vals.SPEED_CFG_CTRL_PERIOD, 4 / 8000, 1e-12, 'nested expression');
  close(vals.MOTOR_CFG_ROTOR_INERTIA, 1.2e-4, 1e-12, 'exponent with suffix');
  assert.ok(unresolved.includes('CURRENT_CFG_FUNC_DUTY_SET'));
  assert.ok(!('MOTOR_SHUNT_TYPE_1_SHUNT' in vals));
});

test('cfg import maps firmware units onto inputs', () => {
  const { vals } = core.parseCfg(fs.readFileSync(fx('sample_cfg.h'), 'utf8'));
  const { inputs, applied } = core.cfgToInputs(vals, core.DEFAULTS);
  assert.equal(inputs.pp, '4');
  assert.equal(inputs.Ld, '4.1');                 // H → mH
  assert.equal(inputs.Jm, '1.2');                 // kg·m² → kg·cm²
  assert.equal(inputs.JL, '0');
  close(+inputs.Irated, 6 * Math.SQRT2, 1e-6, 'Arms → peak');
  assert.equal(inputs.fsw, '8000');               // kHz → Hz
  assert.equal(inputs.fluxMode, 'psia'); assert.equal(inputs.tf, 'pow');
  assert.equal(inputs.curMethod, 'fw'); assert.equal(inputs.fcc, '300'); assert.equal(inputs.zcc, '0.9');
  assert.equal(inputs.ilimC, '150');
  assert.ok(applied.length >= 15);
  const R = core.computeAll(inputs);
  assert.ok(!R.bad);
  close(R.P.psi * Math.sqrt(1.5), 0.0612, 1e-12, 'ψa round trip');
});

test('exported C header uses valid float literals only', () => {
  for (const s of [{}, { piImpl: 'inc', curMethod: 'imc', fcc: '300' }, { tmpTopo: 'tmpL', isnTopo: 'isnH' }]) {
    const R = core.computeAll(s);
    for (const kind of ['c', 'rn']) {
      const lines = core.exportCode(R, kind, '2026-01-01').split('\n').filter(l => l.startsWith('#define') && l.includes('('));
      for (const l of lines) {
        const m = /\(([^()]*)\)/.exec(l); if (!m) continue;
        assert.match(m[1], /^(-?\d+(\.\d+)?(e-?\d+)?f|\d+U?|0x[0-9A-F]+)$/, `bad literal in: ${l}`);
      }
    }
  }
});

test('results match the stored snapshot (catches unintended changes)', () => {
  const snap = JSON.parse(fs.readFileSync(fx('snapshot.json'), 'utf8'));
  for (const [name, { inputs, values }] of Object.entries(snap)) {
    const r = core.computeAll(inputs);
    const now = {
      curKp: r.cur.q.Kp, curKiTs: r.cur.q.Ki * r.P.Ts, curPM: r.cur.q.sw.pm, curBW: r.cur.q.sw.bw, curOS: r.cur.q.met.os,
      spdKp: r.spd.Kp, spdKi: r.spd.Ki, spdPM: r.spd.sw.pm, spdFc: r.spd.sw.fc, spdOS: r.spd.sim.met.os, spdDip: r.spd.sim.dip,
      posPM: r.pos.sw.pm, pllKp: r.obs.pll.Kp, emfKp: r.obs.emf.Kp, loL2: r.obs.lo.l2,
      hwGVA: r.hw.v.GVA, hwVoff: r.hw.v.Voff, hwKv: r.hw.v.Kv, hwTlsb: r.hw.v.tLsb, hwCnt: r.hw.v.cnt };
    for (const [k, v] of Object.entries(values)) close(now[k], v, 1e-6, `${name}.${k}`);
  }
});

test('every input has a default and invalid inputs are reported', () => {
  for (const id of core.FIELD_IDS) assert.equal(typeof core.DEFAULTS[id], 'string', id);
  const R = core.computeAll({ Rs: '0' });
  assert.ok(R.bad && R.bad.includes('Rs'));
});
