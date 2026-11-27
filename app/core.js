/* PMSM tuner — computation core.
 * Pure functions only (no DOM): index.html renders the results, tests/ run the same code under Node.
 * Every input is a string keyed like the form fields; computeAll(inputs) returns the full result object. */
'use strict';
let IN = {};   // inputs of the computation currently running (set by computeAll)

const TAU = 2 * Math.PI;

/* ---------- complex numbers ---------- */
const C = (re, im = 0) => ({ re, im });
const cadd = (a, b) => C(a.re + b.re, a.im + b.im);
const csub = (a, b) => C(a.re - b.re, a.im - b.im);
const cmul = (a, b) => C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return C((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
const cabs = a => Math.hypot(a.re, a.im);
const carg = a => Math.atan2(a.im, a.re);
const cexpj = t => C(Math.cos(t), Math.sin(t));
const cscale = (a, k) => C(a.re * k, a.im * k);
const ONE = C(1, 0);
// zero-order hold factor (1-e^{-jwT})/(jwT)
const zohF = (w, T) => w * T < 1e-9 ? ONE : cdiv(csub(ONE, cexpj(-w * T)), C(0, w * T));
// discrete PI at z=e^{jwT}: (b0 + b1 z^-1)/(1 - z^-1)
const piZ = (w, T, b0, b1) => { const zi = cexpj(-w * T); return cdiv(cadd(C(b0), cscale(zi, b1)), csub(ONE, zi)); };

/* ---------- polynomial roots (Durand–Kerner), coeffs highest first ---------- */
function roots(coef) {
  let c = coef.slice(); while (c.length && Math.abs(c[0]) < 1e-300) c.shift();
  const n = c.length - 1; if (n < 1) return [];
  c = c.map(x => x / c[0]);
  const r = []; let p0 = ONE; for (let i = 0; i < n; i++) { r.push(p0); p0 = cmul(p0, C(0.4, 0.9)); }
  const ev = z => c.reduce((acc, a) => cadd(cmul(acc, z), C(a)), C(0));
  for (let it = 0; it < 500; it++) {
    let mx = 0;
    for (let i = 0; i < n; i++) {
      let den = ONE; for (let j = 0; j < n; j++) if (j !== i) den = cmul(den, csub(r[i], r[j]));
      const d = cdiv(ev(r[i]), den); r[i] = csub(r[i], d); mx = Math.max(mx, cabs(d));
    }
    if (mx < 1e-13) break;
  }
  return r;
}
const polyMul = (a, b) => { const r = new Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => r[i + j] += x * y)); return r; };
const polyAdd = (a, b) => { const n = Math.max(a.length, b.length); const A = Array(n - a.length).fill(0).concat(a), B = Array(n - b.length).fill(0).concat(b); return A.map((x, i) => x + B[i]); };
const eig2 = (a11, a12, a21, a22) => roots([1, -(a11 + a22), a11 * a22 - a12 * a21]);

/* ---------- formatting ---------- */
function fmt(x, p = 4) {
  if (x === undefined || x === null || !isFinite(x)) return '—';
  if (x === 0) return '0';
  const a = Math.abs(x);
  if (a < 1e-3 || a >= 1e6) return x.toExponential(p - 1).replace('e+', 'e');
  return String(+x.toPrecision(p));
}
const cflt = x => { let t = fmt(x, 7); if (!/[.e]/.test(t)) t += '.0'; return t + 'f'; };
const fmtHz = f => !isFinite(f) ? '—' : f >= 1000 ? fmt(f / 1000, 3) + ' kHz' : fmt(f, 3) + ' Hz';
const fmtT = t => !isFinite(t) ? '—' : t < 1e-3 ? fmt(t * 1e6, 3) + ' µs' : t < 1 ? fmt(t * 1e3, 3) + ' ms' : fmt(t, 3) + ' s';
const deg = r => r * 180 / Math.PI;

/* ---------- PI discretization ---------- */
function piCoef(Kp, Ki, T, m) {
  if (m === 'fe') return { b0: Kp, b1: -Kp + Ki * T };
  if (m === 'tu') return { b0: Kp + Ki * T / 2, b1: -Kp + Ki * T / 2 };
  return { b0: Kp + Ki * T, b1: -Kp };
}
const DISC_NAME = { be: '后向欧拉', fe: '前向欧拉', tu: 'Tustin' };

/* ---------- parameters ---------- */
function readP() {
  const v = id => parseFloat(IN[id]);
  const P = {};
  P.type = IN.motorType;
  P.Rs = v('Rs'); P.Ld = v('Ld') * 1e-3; P.Lq = P.type === 'spm' ? P.Ld : v('Lq') * 1e-3;
  P.p = Math.max(1, Math.round(v('pp')));
  const fm = IN.fluxMode, fv = v('fluxVal');
  P.tf = IN.tf; P.kT = P.tf === 'pow' ? Math.sqrt(1.5) : 1;
  P.psi = fm === 'psi' ? fv : fm === 'psia' ? fv / Math.sqrt(1.5) : fm === 'ke' ? fv * Math.SQRT2 / (Math.sqrt(3) * P.p * 1000 * TAU / 60) : fv / (1.5 * P.p);
  P.J = (v('Jm') + v('JL')) * 1e-4; P.Jm = v('Jm') * 1e-4; P.B = v('Bv') || 0;
  P.In = v('Irated'); P.Imax = v('Imax'); P.nN = v('nRated'); P.nMax = v('nMax'); P.Vdc = v('Vdc');
  P.fsw = v('fsw'); P.upd = IN.upd; P.Nd = Math.max(0, Math.min(4, parseFloat(IN.nd) || 0)); P.disc = IN.disc;
  P.piImpl = IN.piImpl === 'rn' ? 'rn' : 'inc'; if (P.piImpl === 'rn') P.disc = 'be';
  P.Ts = P.upd === 'double' ? 1 / (2 * P.fsw) : 1 / P.fsw; P.fs = 1 / P.Ts;
  P.Td = (P.Nd + 0.5) * P.Ts;
  P.Nw = Math.max(1, Math.round(v('Nw'))); P.Tw = P.Nw * P.Ts;
  P.Np = Math.max(1, Math.round(v('Np'))); P.Tp = P.Np * P.Ts;
  P.ff = v('ff') || 0; P.enc = v('enc') || 0; P.spdSrc = IN.spdSrc === 'enc' ? 'enc' : 'pll';
  P.Vmax = P.Vdc / Math.sqrt(3);
  P.Vb = v('Vb') || P.Vmax; P.Ib = v('Ib') || P.Imax; P.nb = v('nb') || P.nMax;
  P.wmb = P.nb * TAU / 60;
  // torque constants
  P.Kt0 = 1.5 * P.p * P.psi;
  const dL = P.Lq - P.Ld;
  if (P.type === 'ipm' && Math.abs(dL) > 1e-12) {
    const Is = P.In;
    let id = (P.psi - Math.sqrt(P.psi * P.psi + 8 * dL * dL * Is * Is)) / (4 * dL);
    if (dL < 0) id = 0; // Ld > Lq: MTPA uses positive id, unusual; fall back to id=0
    const iq = Math.sqrt(Math.max(0, Is * Is - id * id));
    P.idM = id; P.iqM = iq;
    P.TeM = 1.5 * P.p * (P.psi * iq + (P.Ld - P.Lq) * id * iq);
    P.KtM = P.TeM / Is;
  } else { P.idM = 0; P.iqM = P.In; P.TeM = P.Kt0 * P.In; P.KtM = P.Kt0; }
  P.ktMode = P.type === 'ipm' ? IN.ktMode : 'id0';
  P.Kt = P.ktMode === 'mtpa' ? P.KtM : P.Kt0;
  P.TN = P.type === 'ipm' && P.ktMode === 'mtpa' ? P.TeM : P.Kt0 * P.In;
  // firmware-style PI limits, entered in the selected dq convention, stored in amplitude-invariant units
  P.ilimC = (v('ilimC') > 0 ? v('ilimC') : 0.5 * P.Vdc) / P.kT;
  P.ilimS = (v('ilimS') > 0 ? v('ilimS') / P.kT : P.In);
  return P;
}

function validP(P) {
  const bad = [];
  [['Rs', P.Rs], ['Ld', P.Ld], ['Lq', P.Lq], ['ψf', P.psi], ['J', P.J], ['Irated', P.In], ['Imax', P.Imax], ['Udc', P.Vdc], ['fpwm', P.fsw], ['nmax', P.nMax]]
    .forEach(([n, x]) => { if (!(x > 0)) bad.push(n); });
  return bad;
}

/* ---------- frequency sweep analysis ---------- */
function sweep(fmin, fmax, n, Lf, Tf) {
  const out = { f: [], mag: [], ph: [], cl: [] };
  let prev = null, off = 0;
  for (let i = 0; i < n; i++) {
    const f = fmin * Math.pow(fmax / fmin, i / (n - 1)); const w = TAU * f;
    const L = Lf(w); let ph = deg(carg(L));
    if (prev !== null) { while (ph + off - prev > 180) off -= 360; while (ph + off - prev < -180) off += 360; }
    ph += off; prev = ph;
    const T = Tf ? Tf(w, L) : cdiv(L, cadd(ONE, L));
    out.f.push(f); out.mag.push(20 * Math.log10(cabs(L))); out.ph.push(ph); out.cl.push(20 * Math.log10(cabs(T)));
  }
  // margins
  let fc = NaN, pm = NaN, gm = Infinity, fg = NaN, bw = NaN, mr = -Infinity;
  for (let i = 1; i < n; i++) {
    if (isNaN(fc) && out.mag[i - 1] >= 0 && out.mag[i] < 0) {
      const t = out.mag[i - 1] / (out.mag[i - 1] - out.mag[i]);
      fc = out.f[i - 1] * Math.pow(out.f[i] / out.f[i - 1], t);
      const ph = out.ph[i - 1] + t * (out.ph[i] - out.ph[i - 1]);
      pm = 180 + ph; while (pm > 180) pm -= 360; while (pm < -180) pm += 360;
    }
  }
  for (let i = 1; i < n && !isNaN(fc); i++) {
    if (out.f[i] <= fc) continue;
    const a = (out.ph[i - 1] + 180) / 360, b = (out.ph[i] + 180) / 360;
    if (Math.floor(a) !== Math.floor(b)) {
      const lvl = Math.max(Math.floor(a), Math.floor(b)), t = (a - lvl) / (a - b);
      gm = -(out.mag[i - 1] + t * (out.mag[i] - out.mag[i - 1])); fg = out.f[i - 1] * Math.pow(out.f[i] / out.f[i - 1], t); break;
    }
  }
  for (let i = 0; i < n; i++) { mr = Math.max(mr, out.cl[i]); if (isNaN(bw) && out.cl[i] < -3.0103) bw = out.f[i]; }
  Object.assign(out, { fc, pm, gm, fg, bw, mr });
  return out;
}

/* ---------- step metrics ---------- */
function stepMetrics(t, y, ref, t0 = 0, t1 = Infinity) {
  let mx = -Infinity, t10 = NaN, t90 = NaN, ts = NaN;
  const idx = []; for (let i = 0; i < t.length; i++) if (t[i] >= t0 && t[i] <= t1) idx.push(i);
  idx.forEach(i => { mx = Math.max(mx, y[i]); if (isNaN(t10) && y[i] >= 0.1 * ref) t10 = t[i]; if (isNaN(t90) && y[i] >= 0.9 * ref) t90 = t[i]; });
  for (let j = idx.length - 1; j >= 0; j--) { const i = idx[j]; if (Math.abs(y[i] - ref) > 0.02 * Math.abs(ref)) { ts = t[i] - t0; break; } }
  return { os: Math.max(0, (mx - ref) / ref * 100), tr: t90 - t10, ts };
}

/* ============ CURRENT LOOP ============ */
/* Controller design for one axis. L is the inductance the design assumes. */
function currentCtrl(P, L, m, fcc, zc) {
  const wc = TAU * fcc, a = Math.exp(-P.Rs * P.Ts / L), b = (1 - a) / P.Rs;
  let Kp, Ki, b0, b1;
  if (m === 'fw') {
    Kp = 2 * zc * wc * L - P.Rs; Ki = wc * wc * L; ({ b0, b1 } = piCoef(Kp, Ki, P.Ts, P.disc));
  } else if (m === 'zd') {
    const K = 1 - Math.exp(-wc * P.Ts);
    b0 = K / b; b1 = -a * b0;
    // equivalent continuous gains for the chosen discretization
    if (P.disc === 'fe') { Kp = b0; Ki = (b0 + b1) / P.Ts; }
    else if (P.disc === 'tu') { Ki = (b0 + b1) / P.Ts; Kp = b0 - Ki * P.Ts / 2; }
    else { Kp = -b1; Ki = (b0 + b1) / P.Ts; }
  } else { Kp = wc * L; Ki = wc * P.Rs; ({ b0, b1 } = piCoef(Kp, Ki, P.Ts, P.disc)); }
  return { Kp, Ki, b0, b1, Ti: Kp / Ki, KiTI: Ki * P.Ts / Kp };
}
/* Discrete loop of a fixed controller against a plant 1/(Lp·s + Rp) with the configured delay. */
function currentLoop(P, ctrl, Lp, Rp, fref, npts = 500) {
  const Pp = Rp === P.Rs ? P : { ...P, Rs: Rp };
  const pl = plantZ(Pp, Lp);
  // closed-loop characteristic polynomial: den + num, with L(z) = C(z)·G(z)
  const num = polyMul([ctrl.b0, ctrl.b1], pl.num), den = polyMul([1, -1], pl.den);
  const poles = roots(polyAdd(den, num));
  const Lw = w => { const zi = cexpj(-w * P.Ts);
    return cmul(cdiv(cadd(C(ctrl.b0), cscale(zi, ctrl.b1)), csub(ONE, zi)), pl.G(w)); };
  const sw = sweep(Math.max(1, fref / 100), P.fs / 2 * 0.999, npts, Lw);
  return { a: pl.a, b: (1 - pl.a) / Rp, g0: pl.g0, g1: pl.g1, poles, Lw, sw, rmax: Math.max(...poles.map(cabs)) };
}
function designCurrent(P, o = {}) {
  const fcc = o.fcc !== undefined ? o.fcc : parseFloat(IN.fcc), wc = TAU * fcc, m = IN.curMethod, zc = parseFloat(IN.zcc) || 1;
  const ax = L => { const c = currentCtrl(P, L, m, fcc, zc); return { L, ...c, ...currentLoop(P, c, L, P.Rs, fcc, o.npts) }; };
  const d = ax(P.Ld), q = ax(P.Lq);
  if (o.sim !== false) {
    const amp = parseFloat(IN.curStep) || 1;
    d.sim = simCurrent(P, d, amp); q.sim = simCurrent(P, q, amp);
    d.met = stepMetrics(d.sim.t, d.sim.y, amp); q.met = stepMetrics(q.sim.t, q.sim.y, amp);
  }
  return { fcc, wc, m, d, q, fcc45: 1 / (8 * P.Td), fcc60: 1 / (12 * P.Td) };
}

// sampled PI + continuous RL plant; the computed voltage takes effect Nd·Ts after sampling (Nd may be fractional), then ZOH
function simCurrent(P, ax, amp) {
  const Ts = P.Ts, wcl = TAU * Math.max(1, parseFloat(IN.fcc));
  const Tend = Math.max(10 / wcl, 40 * Ts), nS = Math.ceil(Tend / Ts), sub = 20, dt = Ts / sub, dSub = Math.round(P.Nd * sub);
  const ad = Math.exp(-P.Rs * dt / ax.L), bd = (1 - ad) / P.Rs;
  const pi = makePI(P.piImpl, ax.Kp, ax.Ki * Ts, ax.b0, ax.b1, P.piImpl === 'rn' ? P.ilimC : P.Vmax);
  let i = 0, u = 0; const q = [];
  const t = [], y = [];
  for (let k = 0; k < nS; k++) {
    const un = Math.max(-P.Vmax, Math.min(P.Vmax, pi.step(amp - i)));
    q.push([k * sub + dSub, un]);
    for (let s = 0; s < sub; s++) {
      while (q.length && q[0][0] <= k * sub + s) u = q.shift()[1];
      i = ad * i + bd * u; if (s % 4 === 3) { t.push((k * sub + s + 1) * dt); y.push(i); }
    }
  }
  return { t, y, amp };
}

/* PI implementations used by the time-domain simulations
 * inc: u[k] = sat(u[k-1] + b0·e[k] + b1·e[k-1])            (incremental form, output limit only)
 * rn : refi = sat(refi + Ki·T·e); u = sat(Kp·e + refi)     (Renesas motor_pi_ctrl, shared limit, backward Euler) */
function makePI(mode, Kp, KiT, b0, b1, lim) {
  const sat = x => Math.max(-lim, Math.min(lim, x));
  if (mode === 'rn') { let refi = 0; return { step: e => { refi = sat(refi + e * KiT); return sat(Kp * e + refi); } }; }
  let u = 0, ep = 0; return { step: e => { u = sat(u + b0 * e + b1 * ep); ep = e; return u; } };
}

/* ZOH-equivalent RL plant 1/(Ls+R) with an input delay of Nd = m + f samples (modified z-transform)
 * G(z) = z^-m · (g0·z^-1 + g1·z^-2) / (1 − a·z^-1),  g0 = (1−e^{−(1−f)RTs/L})/R,  g1 = (e^{−(1−f)RTs/L} − a)/R
 * Returned as num/den polynomials in z (highest power first) and a frequency-response function. */
function plantZ(P, L) {
  const a = Math.exp(-P.Rs * P.Ts / L), m = Math.floor(P.Nd + 1e-9), f = P.Nd - m;
  const ef = Math.exp(-P.Rs * (1 - f) * P.Ts / L), g0 = (1 - ef) / P.Rs, g1 = f > 1e-9 ? (ef - a) / P.Rs : 0;
  let num, den;
  if (g1 === 0) { num = [g0]; den = [1, -a].concat(Array(m).fill(0)); }          // g0 / (z^m (z − a))
  else { num = [g0, g1]; den = [1, -a].concat(Array(m + 1).fill(0)); }           // (g0 z + g1) / (z^(m+1) (z − a))
  const G = w => { const zi = cexpj(-w * P.Ts); let zm = ONE; for (let k = 0; k < m; k++) zm = cmul(zm, zi);
    return cmul(zm, cdiv(cadd(cscale(zi, g0), cscale(cmul(zi, zi), g1)), csub(ONE, cscale(zi, a)))); };
  return { a, g0, g1, m, f, num, den, G };
}

/* ============ SPEED LOOP ============ */
/* PLL speed estimator running every Ts (forward-Euler form, as used for sensorless angle/speed):
 *   ε = θ − θ̂,  ω̂ = Kp·ε + I,  I += Ki·Ts·ε,  θ̂ += Ts·ω̂
 * Its response from true to estimated speed is ω̂/ω = C(z) / (jω·(1 + C(z)·Ts·z⁻¹/(1 − z⁻¹))), C(z) = Kp + Ki·Ts·z⁻¹/(1 − z⁻¹). */
function pllGains() { const wn = TAU * (parseFloat(IN.fpll) || 25), z = parseFloat(IN.zpll) || 1; return { Kp: 2 * z * wn, Ki: wn * wn, wn, z }; }
function pllResp(P, g, w) {
  const zi = cexpj(-w * P.Ts), integ = cdiv(zi, csub(ONE, zi));
  const Cp = cadd(C(g.Kp), cscale(integ, g.Ki * P.Ts)), A = cscale(integ, P.Ts);
  return cdiv(Cp, cmul(C(0, w), cadd(ONE, cmul(Cp, A))));
}
function speedCtrl(P, cur, o = {}) {
  const m = IN.spdMethod, wf = P.ff > 0 ? TAU * P.ff : Infinity;
  let Kp, Ki, wsc, N, a;
  const Tsig = 1 / cur.wc + (isFinite(wf) ? 1 / wf : 0) + P.Tw + (P.spdSrc === 'pll' ? 2 * pllGains().z / pllGains().wn : 0);
  const fsc = o.fsc !== undefined ? o.fsc : parseFloat(IN.fsc);
  if (m === 'fw') {
    const zs = parseFloat(IN.zsc) || 1;
    wsc = TAU * fsc; Kp = 2 * zs * wsc * P.J / P.Kt; Ki = wsc * wsc * P.J / P.Kt; N = 2 * zs;
  } else if (m === 'so') {
    a = o.aSO !== undefined ? o.aSO : (parseFloat(IN.aSO) || 2.4);
    Kp = P.J / (a * P.Kt * Tsig); const Ti = a * a * Tsig; Ki = Kp / Ti; wsc = 1 / (a * Tsig); N = a;
  } else {
    wsc = TAU * fsc; N = parseFloat(IN.Nz) || 4;
    Kp = P.J * wsc / P.Kt; Ki = Kp * wsc / N;
  }
  const { b0, b1 } = piCoef(Kp, Ki, P.Tw, P.disc);
  return { m, Kp, Ki, b0, b1, wsc, fsc: wsc / TAU, N, a, Tsig, wf };
}
/* Multi-rate frequency response of the speed loop for a fixed controller.
 * plant: { J, Kt, B } may differ from the design values; curQ is the (possibly perturbed) q-axis current loop. */
function speedLoop(P, ctrl, curQ, plant, npts = 500) {
  const wf = ctrl.wf, alpha = isFinite(wf) ? 1 - Math.exp(-wf * P.Tw) : 1;
  // firmware motor_filter_first_order_lpff: y = ka·y[-1] + kb·(x + x[-1]) (Tustin); otherwise y += α(x − y)
  const wT = isFinite(wf) ? wf * P.Tw : Infinity, lpf = P.piImpl === 'rn' && isFinite(wf) ? { type: 'tu', ka: (2 - wT) / (2 + wT), kb: wT / (2 + wT) } : { type: 'iir', alpha };
  const pg = pllGains();
  const Tc = w => { const L = curQ.Lw(w); return cdiv(L, cadd(ONE, L)); };
  const Fz = w => { const zi = cexpj(-w * P.Tw); return lpf.type === 'tu' ? cdiv(cscale(cadd(ONE, zi), lpf.kb), csub(ONE, cscale(zi, lpf.ka))) : cdiv(C(alpha), csub(ONE, cscale(zi, 1 - alpha))); };
  const Mech = w => cdiv(C(plant.Kt), C(plant.B, plant.J * w));
  const fwd = w => cmul(cmul(cmul(piZ(w, P.Tw, ctrl.b0, ctrl.b1), zohF(w, P.Tw)), Tc(w)), Mech(w));
  // speed measurement: PLL estimate (sensorless) or position difference Δθ/Tω (encoder), then the low-pass filter
  const meas = w => P.spdSrc === 'pll' ? pllResp(P, pg, w) : zohF(w, P.Tw);
  const fb = w => cmul(meas(w), Fz(w));
  const Lw = w => cmul(fwd(w), fb(w));
  const Tw = w => cdiv(fwd(w), cadd(ONE, Lw(w)));
  const sw = sweep(Math.max(0.1, ctrl.fsc / 100), Math.min(P.fs / 2, 1 / (2 * P.Tw)) * 0.999, npts, Lw, w => Tw(w));
  return { alpha, lpf, pll: pg, Lw, Tw, sw };
}
function designSpeed(P, cur, o = {}) {
  const c = speedCtrl(P, cur, o);
  const R = { ...c, ...speedLoop(P, c, cur.q, P, o.npts) };
  if (o.sim !== false) R.sim = simCascade(P, cur, R, null, 'speed');
  return R;
}

/* ============ POSITION LOOP ============ */
function designPos(P, cur, spd) {
  const fpc = parseFloat(IN.fpc), Kpp = TAU * fpc;
  const vff = (parseFloat(IN.vff) || 0) / 100;
  const Kaff = P.J / P.Kt;
  const Lw = w => cmul(cmul(C(Kpp), zohF(w, P.Tp)), cmul(spd.Tw(w), cdiv(ONE, C(0, w))));
  const sw = sweep(Math.max(0.05, fpc / 100), Math.min(1 / (2 * P.Tp), spd.fsc * 20), 500, Lw);
  const R = { fpc, Kpp, vff, Kaff, Lw, sw };
  R.sim = simCascade(P, cur, spd, R, 'pos');
  return R;
}

/* ============ AUTO-TUNE, TRADE-OFF, ROBUSTNESS ============ */
/* All three use the frequency-domain analysis only (no time simulation), so they stay fast. */
const loopMetrics = (loop, R) => loop === 'cur'
  ? { pm: Math.min(R.d.sw.pm, R.q.sw.pm), gm: Math.min(R.d.sw.gm, R.q.sw.gm), bw: R.q.sw.bw, fc: R.q.sw.fc, rmax: Math.max(R.d.rmax, R.q.rmax) }
  : { pm: R.sw.pm, gm: R.sw.gm, bw: R.sw.bw, fc: R.sw.fc, rmax: 0 };
function tuneSpace(loop, P, cur) {
  if (loop === 'cur') return { key: 'fcc', lo: Math.max(1, P.fs / 5000), hi: P.fs / 3, larger: true };
  if (IN.spdMethod === 'so') return { key: 'aSO', lo: 1.2, hi: 8, larger: false };
  return { key: 'fsc', lo: Math.max(0.05, cur.fcc / 2000), hi: Math.min(cur.fcc, 1 / (2 * P.Tw)), larger: true };
}
function evalTune(loop, P, cur, key, x) {
  const R = loop === 'cur' ? designCurrent(P, { fcc: x, sim: false, npts: 300 }) : designSpeed(P, cur, { [key]: x, sim: false, npts: 300 });
  return loopMetrics(loop, R);
}
/* plant variation ranges from the inputs; each list contains 1 (nominal) */
function rbRanges() {
  const g = id => parseFloat(IN[id]), uniq = a => a.filter((x, i) => a.indexOf(x) === i);
  const lv = (g('rbL') || 0) / 100, kv = (g('rbKt') || 0) / 100;
  return { Ls: lv > 0 ? [1 - lv, 1, 1 + lv] : [1], Rs: uniq([g('rbRlo') || 1, 1, g('rbRhi') || 1]),
    Js: uniq([g('rbJlo') || 1, 1, g('rbJhi') || 1]), Ks: kv > 0 ? [1 - kv, 1, 1 + kv] : [1] };
}
/* worst margins over all plant corners for the controller that the tuning parameter x produces */
function worstAt(loop, P, cur, key, x, rr) {
  let pm = Infinity, gm = Infinity, rmax = 0, bw = NaN;
  if (loop === 'cur') {
    const c = designCurrent(P, { fcc: x, sim: false, npts: 200 });
    for (const ls of rr.Ls) for (const rs of rr.Rs) for (const ax of ['d', 'q']) {
      const l = currentLoop(P, c[ax], (ax === 'd' ? P.Ld : P.Lq) * ls, P.Rs * rs, x, 200);
      pm = Math.min(pm, l.sw.pm); if (isFinite(l.sw.gm)) gm = Math.min(gm, l.sw.gm); rmax = Math.max(rmax, l.rmax);
    }
    bw = c.q.sw.bw;
  } else {
    const sc = speedCtrl(P, cur, { [key]: x });
    for (const ls of rr.Ls) {
      const q = currentLoop(P, cur.q, P.Lq * ls, P.Rs, cur.fcc, 150);
      for (const js of rr.Js) for (const ks of rr.Ks) {
        const l = speedLoop(P, sc, q, { J: P.J * js, Kt: P.Kt * ks, B: P.B }, 200);
        pm = Math.min(pm, l.sw.pm); if (isFinite(l.sw.gm)) gm = Math.min(gm, l.sw.gm);
        if (ls === 1 && js === 1 && ks === 1) bw = l.sw.bw;
      }
    }
  }
  return { pm, gm, rmax, bw };
}
/* Bandwidth–margin trade-off: margins over the tuning parameter, log-spaced. */
function tradeoff(inputs, loop, n = 40) {
  IN = { ...DEFAULTS, ...inputs };
  const P = readP(); if (validP(P).length) return null;
  const cur = designCurrent(P, { sim: false, npts: 300 });
  const sp = tuneSpace(loop, P, cur), xs = [], pts = [];
  for (let i = 0; i < n; i++) {
    const x = sp.lo * Math.pow(sp.hi / sp.lo, i / (n - 1));
    xs.push(x); pts.push(evalTune(loop, P, cur, sp.key, x));
  }
  return { ...sp, xs, pm: pts.map(p => p.pm), gm: pts.map(p => p.gm), bw: pts.map(p => p.bw) };
}
/* Largest bandwidth (or smallest a for symmetric optimum) that still meets the phase/gain margin targets.
 * mode 'worst' requires the targets at every plant corner of the robustness ranges, not only at nominal values. */
function autoTune(inputs, loop, pmTarget, gmTarget, mode) {
  IN = { ...DEFAULTS, ...inputs };
  const P = readP(); if (validP(P).length) return null;
  const cur = designCurrent(P, { sim: false, npts: 300 });
  const sp = tuneSpace(loop, P, cur), worst = mode === 'worst', rr = rbRanges();
  const metr = x => worst ? worstAt(loop, P, cur, sp.key, x, rr) : evalTune(loop, P, cur, sp.key, x);
  const ok = x => { const m = metr(x); return m.pm >= pmTarget && (!isFinite(m.gm) || m.gm >= gmTarget) && m.rmax < 1 && isFinite(m.pm); };
  const n = worst ? 32 : 48, xs = []; for (let i = 0; i < n; i++) xs.push(sp.lo * Math.pow(sp.hi / sp.lo, i / (n - 1)));
  // scan from the conservative end; the first infeasible point after a feasible run brackets the limit
  const order = sp.larger ? xs : xs.slice().reverse();
  let good = null, bad = null;
  for (const x of order) { if (ok(x)) good = x; else if (good !== null) { bad = x; break; } }
  if (good === null) return { key: sp.key, value: null };
  if (bad !== null) for (let i = 0; i < (worst ? 18 : 30); i++) { const mid = Math.sqrt(good * bad); if (ok(mid)) good = mid; else bad = mid; }
  const m = evalTune(loop, P, cur, sp.key, good), w = worst ? worstAt(loop, P, cur, sp.key, good, rr) : null;
  return { key: sp.key, value: good, mode: worst ? 'worst' : 'nom', ...m, worstPm: w ? w.pm : undefined, worstGm: w ? w.gm : undefined };
}
/* Fixed controller (designed with the nominal parameters) against a plant whose parameters vary over the given ranges. */
function robustness(inputs, loop) {
  IN = { ...DEFAULTS, ...inputs };
  const P = readP(); if (validP(P).length) return null;
  const { Ls, Rs, Js, Ks } = rbRanges();
  const cur = designCurrent(P, { sim: false });
  const rows = [];
  if (loop === 'cur') {
    for (const ls of Ls) for (const rs of Rs) {
      const d = currentLoop(P, cur.d, P.Ld * ls, P.Rs * rs, cur.fcc, 300), q = currentLoop(P, cur.q, P.Lq * ls, P.Rs * rs, cur.fcc, 300);
      rows.push({ L: ls, R: rs, pm: Math.min(d.sw.pm, q.sw.pm), gm: Math.min(d.sw.gm, q.sw.gm), bw: q.sw.bw, rmax: Math.max(d.rmax, q.rmax) });
    }
  } else {
    const spd = speedCtrl(P, cur);
    for (const ls of Ls) {
      const q = currentLoop(P, cur.q, P.Lq * ls, P.Rs, cur.fcc, 200);
      for (const js of Js) for (const ks of Ks) {
        const l = speedLoop(P, spd, q, { J: P.J * js, Kt: P.Kt * ks, B: P.B }, 300);
        rows.push({ L: ls, J: js, Kt: ks, pm: l.sw.pm, gm: l.sw.gm, bw: l.sw.bw, fc: l.sw.fc });
      }
    }
  }
  const fin = a => a.filter(isFinite);
  const worst = { pm: Math.min(...rows.map(r => r.pm)), gm: Math.min(...fin(rows.map(r => r.gm)).concat([Infinity])), bwLo: Math.min(...rows.map(r => r.bw)), bwHi: Math.max(...rows.map(r => r.bw)) };
  return { rows, worst, nominal: rows.find(r => r.L === 1 && (r.R === 1 || r.R === undefined) && (r.J === 1 || r.J === undefined) && (r.Kt === 1 || r.Kt === undefined)) };
}

/* multi-rate cascade simulation: current (Ts) / speed (Tw) / position (Tp) */
function simCascade(P, cur, spd, pos, mode) {
  const Ts = P.Ts, q = cur.q;
  const zoz = Math.max(1e-3, spd.Ki / spd.Kp);
  let Tend, Tload = Infinity, Tl = 0;
  const step = mode === 'speed' ? (parseFloat(IN.spdStep) || 100) * TAU / 60 : (parseFloat(IN.posStep) || 10) * Math.PI / 180;
  if (mode === 'speed') {
    const Tst = Math.min(4, 6 / spd.wsc + 3 / zoz);
    const lp = (parseFloat(IN.loadPct) || 0) / 100;
    if (lp > 0) { Tload = Tst; Tl = lp * P.TN; Tend = 2 * Tst; } else Tend = Tst;
  } else Tend = Math.min(6, 10 / pos.Kpp + 3 / zoz);
  let sub = 8; let nS = Math.ceil(Tend / Ts);
  if (nS * sub > 600000) sub = Math.max(2, Math.floor(600000 / nS));
  if (nS > 600000) nS = 600000;
  const dt = Ts / sub;
  const ad = Math.exp(-P.Rs * dt / q.L), bd = (1 - ad) / P.Rs, dSub = Math.round(P.Nd * sub);
  const rn = P.piImpl === 'rn';
  const piW = makePI(P.piImpl, spd.Kp, spd.Ki * P.Tw, spd.b0, spd.b1, rn ? P.ilimS : P.Imax);
  const piC = makePI(P.piImpl, q.Kp, q.Ki * Ts, q.b0, q.b1, rn ? P.ilimC : P.Vmax * 2);
  let iq = 0, w = 0, th = 0, thPrevW = 0, wf = 0, wmPrev = 0, iqRef = 0, wRef = 0;
  const pll = spd.pll || pllGains(); let thHat = 0, pllI = 0, wHat = 0;   // PLL on the electrical angle, runs every Ts
  const dq = [];
  let u = 0;
  const t = [], y = [], y2 = [], r = [];
  const dec = Math.max(1, Math.floor(nS * sub / 1500));
  let cnt = 0;
  for (let k = 0; k < nS; k++) {
    const tk = k * Ts;
    // position loop
    if (mode === 'pos' && k % P.Np === 0) {
      wRef = pos.Kpp * (step - th); const wl = P.nMax * TAU / 60; wRef = Math.max(-wl, Math.min(wl, wRef));
    }
    if (mode === 'speed') wRef = step;
    // speed estimate: PLL (sensorless) every Ts, or position difference at the speed-loop rate
    if (P.spdSrc === 'pll') { const eps = P.p * th - thHat; wHat = pll.Kp * eps + pllI; pllI += pll.Ki * Ts * eps; thHat += Ts * wHat; }
    // speed loop
    if (k % P.Nw === 0) {
      const wm = P.spdSrc === 'pll' ? wHat / P.p : (th - thPrevW) / P.Tw; thPrevW = th;
      if (spd.lpf.type === 'tu') { wf = spd.lpf.ka * wf + spd.lpf.kb * (wm + wmPrev); wmPrev = wm; } else wf += spd.alpha * (wm - wf);
      iqRef = piW.step(wRef - wf);
    }
    // current loop (q axis) with decoupling feed-forward ωe·ψf
    const ffv = P.p * w * P.psi;
    const uo = Math.max(-P.Vmax, Math.min(P.Vmax, piC.step(iqRef - iq) + ffv));
    dq.push([k * sub + dSub, uo]);
    for (let s = 0; s < sub; s++) {
      while (dq.length && dq[0][0] <= k * sub + s) u = dq.shift()[1];
      const emf = P.p * w * P.psi;
      iq = ad * iq + bd * (u - emf);
      const TL = (tk + s * dt) >= Tload ? Tl * Math.sign(step || 1) : 0;
      const Te = P.Kt * iq;
      w += dt * (Te - TL - P.B * w) / P.J;
      th += dt * w;
      if (++cnt % dec === 0) {
        t.push(tk + (s + 1) * dt);
        if (mode === 'speed') { y.push(w * 60 / TAU); r.push(wRef * 60 / TAU); }
        else { y.push(th * 180 / Math.PI); r.push(step * 180 / Math.PI); }
        y2.push(iq);
      }
    }
  }
  const ref = mode === 'speed' ? step * 60 / TAU : step * 180 / Math.PI;
  const met = stepMetrics(t, y, ref, 0, Tload);
  let dip = NaN, rec = NaN;
  if (mode === 'speed' && isFinite(Tload)) {
    let mn = Infinity, ti = 0;
    for (let i = 0; i < t.length; i++) if (t[i] >= Tload) { if (y[i] < mn) { mn = y[i]; ti = i; } }
    dip = ref - mn;
    for (let i = t.length - 1; i >= 0; i--) { if (t[i] < Tload) break; if (Math.abs(y[i] - ref) > 0.02 * Math.abs(ref)) { rec = t[i] - Tload; break; } }
  }
  return { t, y, y2, r, ref, met, dip, rec, Tload, Tl };
}

/* ============ OBSERVERS ============ */
function designObs(P) {
  const R = {};
  const wn = TAU * parseFloat(IN.fpll), z = parseFloat(IN.zpll) || 0.707;
  R.pll = { wn, z, Kp: 2 * z * wn, Ki: wn * wn };
  R.pll.bw = wn * Math.sqrt(1 + 2 * z * z + Math.sqrt(Math.pow(1 + 2 * z * z, 2) + 1)) / TAU;
  R.pll.poles = eig2(1 - P.Ts * R.pll.Kp, -P.Ts, R.pll.Ki * P.Ts, 1);
  R.pll.Emin = P.p * (P.nN * 0.05 * TAU / 60) * P.psi;

  const wo = TAU * parseFloat(IN.femf), zo = parseFloat(IN.zemf) || 0.707, L = P.Lq;
  const e = { wo, z: zo, Kp: 2 * zo * wo * L - P.Rs, Ki: wo * wo * L };
  const weN = P.p * P.nN * TAU / 60;
  const H = cdiv(C(e.Ki, e.Kp * weN), C(e.Ki - L * weN * weN, (P.Rs + e.Kp) * weN));
  e.gainN = cabs(H); e.phN = deg(carg(H)); e.weN = weN;
  e.poles = eig2(1 - P.Ts * (P.Rs + e.Kp) / L, -P.Ts / L, e.Ki * P.Ts, 1);
  R.emf = e;

  const eta = parseFloat(IN.smoEta) || 1.5, flp = parseFloat(IN.smoLpf) || 1000;
  const weMax = P.p * P.nMax * TAU / 60;
  const s = { eta, k: eta * weMax * P.psi, flp, wlp: TAU * flp, weMax };
  s.lagN = deg(Math.atan(weN / s.wlp)); s.lagMax = deg(Math.atan(weMax / s.wlp));
  s.ampN = 1 / Math.sqrt(1 + Math.pow(weN / s.wlp, 2));
  s.alpha = 1 - Math.exp(-s.wlp * P.Ts);
  s.ripple = s.k * P.Ts / P.Ld;
  R.smo = s;

  const wl = TAU * parseFloat(IN.flo), zl = parseFloat(IN.zlo) || 0.707;
  const lo = { wo: wl, z: zl, l1: 2 * zl * wl - P.B / P.J, l2: P.J * wl * wl };
  lo.poles = eig2(1 - P.Tw * (P.B / P.J + lo.l1), -P.Tw / P.J, P.Tw * lo.l2, 1);
  lo.Kcomp = 1 / P.Kt;
  R.lo = lo;
  return R;
}

/* ============ HARDWARE ============ */
/* SVG drawing helpers for the basic sampling circuits */
const SV = (() => {
  const t = (x, y, s, c = 'lb', a = 'middle') => `<text x="${x}" y="${y}" class="${c}" text-anchor="${a}">${s}</text>`;
  const w = (...p) => `<polyline class="w" points="${p.map(q => q.join(',')).join(' ')}"/>`;
  const dot = (x, y) => `<circle class="dt" cx="${x}" cy="${y}" r="2.8"/>`;
  function R(x1, y1, x2, y2, lab, val, o = {}) {
    const vert = x1 === x2, mx = (x1 + x2) / 2, my = (y1 + y2) / 2, L = 30, W = 11, d = Math.sign(vert ? y2 - y1 : x2 - x1) || 1;
    let s = '';
    if (vert) {
      s += w([x1, y1], [x1, my - d * L / 2]) + w([x1, my + d * L / 2], [x2, y2]) + `<rect class="cb" x="${mx - W / 2}" y="${my - L / 2}" width="${W}" height="${L}"/>`;
      const sd = o.left ? -1 : 1, a = o.left ? 'end' : 'start';
      s += t(mx + sd * 11, my - 2, lab, 'lb', a) + t(mx + sd * 11, my + 11, val, 'lv', a);
    } else {
      s += w([x1, y1], [mx - d * L / 2, y1]) + w([mx + d * L / 2, y1], [x2, y2]) + `<rect class="cb" x="${mx - L / 2}" y="${my - W / 2}" width="${L}" height="${W}"/>`;
      if (o.above) s += t(mx, my - 22, lab, 'lb') + t(mx, my - 10, val, 'lv');
      else if (o.below) s += t(mx, my + 21, lab, 'lb') + t(mx, my + 33, val, 'lv');
      else s += t(mx, my - 10, lab, 'lb') + t(mx, my + 20, val, 'lv');
    }
    if (o.ntc) s += `<path class="w" d="M${mx - 11},${my + 14} L${mx + 11},${my - 14} h6"/>` + '';
    return s;
  }
  function Cp(x1, y1, x2, y2, lab, val, o = {}) {
    const vert = x1 === x2, mx = (x1 + x2) / 2, my = (y1 + y2) / 2, g = 4, P = 18;
    if (vert) return w([x1, y1], [x1, my - g]) + w([x1, my + g], [x2, y2]) + `<path class="cpl" d="M${mx - P / 2},${my - g}h${P}M${mx - P / 2},${my + g}h${P}"/>`
      + t(mx + (o.left ? -13 : 13), my - 1, lab, 'lb', o.left ? 'end' : 'start') + t(mx + (o.left ? -13 : 13), my + 12, val, 'lv', o.left ? 'end' : 'start');
    return w([x1, y1], [mx - g, y1]) + w([mx + g, y1], [x2, y2]) + `<path class="cpl" d="M${mx - g},${my - P / 2}v${P}M${mx + g},${my - P / 2}v${P}"/>`
      + (o.below ? t(mx, my + 24, lab) + t(mx, my + 36, val, 'lv') : t(mx, my - 13, lab) + t(mx, my + 23, val, 'lv'));
  }
  const gnd = (x, y) => w([x, y], [x, y + 8]) + `<path class="w" d="M${x - 9},${y + 8}h18M${x - 6},${y + 12}h12M${x - 3},${y + 16}h6"/>`;
  const vs = (x, y, s) => `<path class="w" d="M${x - 9},${y}h18"/>` + t(x, y - 6, s);
  function op(x, y, lab, plusTop) {
    const a = plusTop ? '+' : '−', b = plusTop ? '−' : '+';
    return `<path class="cb" d="M${x},${y - 38} L${x},${y + 38} L${x + 66},${y} Z"/>` + t(x + 8, y - 16, a, 'pm', 'start') + t(x + 8, y + 24, b, 'pm', 'start') + t(x + 24, y + 4, lab, 'lv', 'start');
  }
  const port = (x, y, s, a = 'start') => `<circle class="pt" cx="${x}" cy="${y}" r="4"/>` + t(x + (a === 'start' ? 9 : -9), y + 4, s, 'lb', a);
  const box = (x, y, bw, bh, l1, l2 = '') => `<rect class="cb" x="${x}" y="${y}" width="${bw}" height="${bh}" rx="4"/>` + t(x + bw / 2, y + bh / 2 - (l2 ? 3 : -4), l1) + (l2 ? t(x + bw / 2, y + bh / 2 + 13, l2, 'lv') : '');
  const svg = (W, H, body) => `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="基本单元电路">${body}</svg>`;
  return { t, w, dot, R, Cp, gnd, vs, op, port, box, svg };
})();

/* field: [key, label, unit, default, boardRef] — ids are `${topoId}_${key}` */
const UNIT = { 'mΩ': 1e-3, 'Ω': 1, 'kΩ': 1e3, 'MΩ': 1e6, 'pF': 1e-12, 'nF': 1e-9, 'µF': 1e-6, 'mV/A': 1e-3, 'V': 1, '': 1, 'K': 1, '°C': 1, '°C/V': 1, 'mV/°C': 1e-3, 'bit': 1, 'LSB': 1, 'V/V': 1 };
const DIFF_FIELDS = [['R1', '采样电阻', 'mΩ', 20, 'R83 (V相 R840)'], ['R2', '同相侧串联电阻', 'Ω', 10, 'R84'], ['R3', '反相侧串联电阻', 'Ω', 10, 'R85'],
  ['R4', '反相输入电阻 Rg−', 'kΩ', 2.2, 'R89'], ['R5', '同相输入电阻 Rg+', 'kΩ', 2.2, 'R90'], ['R6', '反馈电阻 Rf', 'kΩ', 10, 'R92'],
  ['R7', '偏置上拉 (接 Vs)', 'kΩ', 20, 'R111'], ['R8', '偏置下拉 (接 AGND)', 'kΩ', 20, 'R91'], ['Vs', '偏置电源', 'V', 5, 'AP5V'],
  ['C2', '反馈电容', 'pF', 10, 'C92'], ['R9', '输出 RC 电阻', 'Ω', 470, 'R93'], ['C3', '输出 RC 电容', 'pF', 220, 'C93']];
const diffDiagram = (v, top) => SV.svg(620, 340,
  SV.t(50, 70, top) + SV.w([50, 78], [50, 155]) + SV.R(50, 155, 50, 195, 'R1', v.R1, { left: true }) + SV.w([50, 195], [50, 225]) + SV.gnd(50, 225) + SV.dot(50, 155) + SV.dot(50, 195) +
  SV.w([50, 155], [70, 155]) + SV.R(70, 155, 130, 155, 'R2', v.R2, { above: true }) + SV.w([130, 155], [160, 155]) + SV.R(160, 155, 220, 155, 'R5', v.R5, { above: true }) + SV.w([220, 155], [330, 155]) +
  SV.w([50, 195], [70, 195]) + SV.R(70, 195, 130, 195, 'R3', v.R3, { below: true }) + SV.w([130, 195], [160, 195]) + SV.R(160, 195, 220, 195, 'R4', v.R4, { below: true }) + SV.w([220, 195], [330, 195]) +
  SV.dot(270, 155) + SV.w([270, 155], [270, 100]) + SV.dot(270, 100) + SV.R(270, 100, 270, 38, 'R7', v.R7) + SV.vs(270, 38, 'Vs ' + v.Vs) + SV.R(270, 100, 200, 100, 'R8', v.R8) + SV.gnd(200, 100) +
  SV.dot(300, 195) + SV.w([300, 195], [300, 255]) + SV.dot(300, 255) + SV.R(300, 255, 430, 255, 'R6', v.R6, { above: true }) + SV.w([430, 255], [430, 175]) +
  SV.w([300, 255], [300, 285]) + SV.Cp(300, 285, 430, 285, 'C2', v.C2, { below: true }) + SV.w([430, 285], [430, 255]) + SV.dot(430, 255) +
  SV.op(330, 175, 'U1', true) + SV.w([396, 175], [440, 175]) + SV.dot(430, 175) + SV.R(440, 175, 500, 175, 'R9', v.R9) + SV.w([500, 175], [556, 175]) +
  SV.dot(520, 175) + SV.Cp(520, 175, 520, 235, 'C3', v.C3) + SV.gnd(520, 235) + SV.port(560, 175, 'ADC') +
  SV.t(145, 151, 'V+', 'lv') + SV.t(145, 191, 'V−', 'lv'));
// general op-amp: Vout = a·V+ − b·V− + c  (inputs through Rp into + node biased by Thevenin (Vb, Rb), Rn into − node, Rf feedback)
function opLin(Rp, Rbth, Vbth, Rn, Rf) {
  const k = isFinite(Rbth) ? Rbth / (Rp + Rbth) : 1, kb = isFinite(Rbth) ? Rp / (Rp + Rbth) : 0;
  return { a: (1 + Rf / Rn) * k, b: Rf / Rn, c: (1 + Rf / Rn) * kb * Vbth };
}
const par = (a, b) => a * b / (a + b);
const rcRows = (r, c) => {
  if (!(r > 0 && c > 0)) return [];
  return [['输出 RC 截止', '1/(2π·R·C)', 1 / (TAU * r * c) / 1e3, 'kHz'], ['输出 RC 时间常数', 'R·C', r * c * 1e6, 'µs']];
};

const HWG = [
  { id: 'adc', title: 'MCU ADC', sub: 'RX24T S12AD：12 位，AVCC0 / VREFH0 = 5 V', common: [
    ['adcBits', 'ADC 分辨率 N', 'bit', 12], ['adcVref', 'ADC 参考电压 Vref', 'V', 5], ['adcClk', 'ADC 时钟 fAD', '', 40, 'MHz'], ['tAd', '电流采样窗口 tAD', '', 6, 'µs']] },
  { id: 'tol', title: '容差与运放参数', sub: '用于下面各采样链路的最坏情况分析（所有容差角点组合）', common: [
    ['tolR', '普通电阻精度', '%', 1], ['tolSh', '采样电阻精度', '%', 1], ['tolV', '基准 / 偏置电源精度', '%', 1],
    ['tolSc', 'IPM CIN 阈值精度', '%', 5], ['tolNtc', 'NTC R25 精度', '%', 1], ['tolB', 'NTC B 值精度', '%', 1],
    ['opVos', '运放失调电压 (最大)', 'mV', 2], ['opSw', '运放输出摆幅余量 (距电源轨)', 'mV', 50], ['opGbw', '运放增益带宽积', 'MHz', 10]] },
  { id: 'isn', title: '相电流采样', sub: '输出 Vadc = Voff + G·i，固件 i = (OFFSET − ad)·LSB/G', topos: [
    { id: 'isnA', name: '分流电阻 + 差分放大 + 分压偏置 (本板)', board: 'U10-A，U 相；V 相同构', fields: DIFF_FIELDS,
      draw: v => diffDiagram(v, '下桥 N 端 (NU)'),
      calc: g => { const o = opLin(g.R5 + g.R2, par(g.R7, g.R8), g.Vs * g.R8 / (g.R7 + g.R8), g.R4 + g.R3, g.R6);
        return { Gv: o.a, cm: o.a - o.b, Voff: o.c, ng: 1 + g.R6 / (g.R4 + g.R3), rsh: g.R1, rc: [g.R9, g.C3], fb: [g.R6, g.C2], rows: [
          ['同相支路系数 a', '(1+R6/(R4+R3))·Rth/(R5+R2+Rth)，Rth=R7∥R8', o.a, 'V/V'],
          ['反相支路系数 b', 'R6/(R4+R3)', o.b, 'V/V'],
          ['共模增益 a−b', '≈0 表示电阻匹配', o.a - o.b, 'V/V', Math.abs(o.a - o.b) < 0.01 * o.a ? 'ok' : 'warn'],
          ['偏置等效电阻匹配', 'R7∥R8 = R6 ?', par(g.R7, g.R8) / 1e3, 'kΩ', Math.abs(par(g.R7, g.R8) / g.R6 - 1) < 0.02 ? 'ok' : 'warn']] }; } },
    { id: 'isnS', name: '单电阻母线采样 + 差分放大 (1-shunt)', board: '', fields: DIFF_FIELDS.map(f => [...f.slice(0, 4), '']),
      draw: v => diffDiagram(v, '逆变器 N 母线'),
      calc: g => { const o = opLin(g.R5 + g.R2, par(g.R7, g.R8), g.Vs * g.R8 / (g.R7 + g.R8), g.R4 + g.R3, g.R6);
        return { Gv: o.a, cm: o.a - o.b, Voff: o.c, ng: 1 + g.R6 / (g.R4 + g.R3), rsh: g.R1, rc: [g.R9, g.C3], fb: [g.R6, g.C2], note: '单电阻需在一个 PWM 周期内两次采样重构三相电流，采样窗口要求更严格。', rows: [
          ['同相支路系数 a', '(1+R6/(R4+R3))·Rth/(R5+R2+Rth)', o.a, 'V/V'], ['共模增益 a−b', '≈0 表示电阻匹配', o.a - o.b, 'V/V']] }; } },
    { id: 'isnI', name: '分流电阻 + 反相放大 + 同相端偏置', board: '', fields: [
      ['R1', '采样电阻', 'mΩ', 20], ['R2', '输入电阻 Rg', 'kΩ', 2.2], ['R3', '反馈电阻 Rf', 'kΩ', 10], ['R4', '偏置上拉 (接 Vs)', 'kΩ', 10], ['R5', '偏置下拉', 'kΩ', 10], ['Vs', '偏置电源', 'V', 5], ['R6', '输出 RC 电阻', 'Ω', 470], ['C1', '输出 RC 电容', 'pF', 220]],
      draw: v => SV.svg(620, 290,
        SV.t(50, 70, '下桥 N 端') + SV.w([50, 78], [50, 135]) + SV.R(50, 135, 50, 195, 'R1', v.R1, { left: true }) + SV.gnd(50, 195) + SV.dot(50, 135) + SV.w([50, 135], [110, 135]) + SV.R(110, 135, 200, 135, 'R2', v.R2) + SV.w([200, 135], [330, 135]) + SV.dot(280, 135) +
        SV.w([280, 135], [280, 70]) + SV.R(280, 70, 450, 70, 'R3', v.R3) + SV.w([450, 70], [450, 155]) +
        SV.w([330, 175], [230, 175]) + SV.dot(230, 175) + SV.R(230, 175, 230, 250, 'R5', v.R5, { left: true }) + SV.gnd(230, 250) + SV.R(230, 175, 140, 175, 'R4', v.R4) + SV.port(135, 175, 'Vs ' + v.Vs, 'end') +
        SV.op(330, 155, 'U1', false) + SV.w([396, 155], [470, 155]) + SV.dot(450, 155) + SV.R(470, 155, 530, 155, 'R6', v.R6) + SV.w([530, 155], [556, 155]) + SV.dot(545, 155) + SV.Cp(545, 155, 545, 215, 'C1', v.C1) + SV.gnd(545, 215) + SV.port(565, 155, 'ADC')),
      calc: g => { const o = opLin(0, Infinity, 0, g.R2, g.R3); const Vp = g.Vs * g.R5 / (g.R4 + g.R5);
        return { Gv: -g.R3 / g.R2, Voff: Vp * (1 + g.R3 / g.R2), ng: 1 + g.R3 / g.R2, rsh: g.R1, rc: [g.R6, g.C1], rows: [
          ['同相端偏置', 'Vp = Vs·R5/(R4+R5)', Vp, 'V'], ['电压增益', '−R3/R2', -g.R3 / g.R2, 'V/V'], ['输出零点', 'Vp·(1+R3/R2)', Vp * (1 + g.R3 / g.R2), 'V']] }; } },
    { id: 'isnN', name: '分流电阻 + 同相放大 + 电阻抬升偏置', board: '', fields: [
      ['R1', '采样电阻', 'mΩ', 20], ['R2', '采样串联电阻', 'kΩ', 1], ['R3', '抬升电阻 (接 Vs)', 'kΩ', 10], ['Vs', '抬升电源', 'V', 5], ['R4', '增益电阻 (接地)', 'kΩ', 1], ['R5', '反馈电阻', 'kΩ', 4], ['R6', '输出 RC 电阻', 'Ω', 470], ['C1', '输出 RC 电容', 'pF', 220]],
      draw: v => SV.svg(620, 300,
        SV.t(50, 70, '下桥 N 端') + SV.w([50, 78], [50, 135]) + SV.R(50, 135, 50, 195, 'R1', v.R1, { left: true }) + SV.gnd(50, 195) + SV.dot(50, 135) + SV.w([50, 135], [110, 135]) + SV.R(110, 135, 200, 135, 'R2', v.R2) + SV.w([200, 135], [330, 135]) + SV.dot(270, 135) + SV.R(270, 135, 270, 50, 'R3', v.R3) + SV.vs(270, 50, 'Vs ' + v.Vs) +
        SV.w([330, 175], [240, 175]) + SV.dot(240, 175) + SV.R(240, 175, 240, 255, 'R4', v.R4, { left: true }) + SV.gnd(240, 255) +
        SV.dot(290, 175) + SV.w([290, 175], [290, 225]) + SV.R(290, 225, 450, 225, 'R5', v.R5) + SV.w([450, 225], [450, 155]) +
        SV.op(330, 155, 'U1', true) + SV.w([396, 155], [470, 155]) + SV.dot(450, 155) + SV.R(470, 155, 530, 155, 'R6', v.R6) + SV.w([530, 155], [556, 155]) + SV.dot(545, 155) + SV.Cp(545, 155, 545, 215, 'C1', v.C1) + SV.gnd(545, 215) + SV.port(565, 155, 'ADC')),
      calc: g => { const k = g.R3 / (g.R2 + g.R3), A = 1 + g.R5 / g.R4;
        return { Gv: k * A, Voff: g.Vs * g.R2 / (g.R2 + g.R3) * A, ng: A, rsh: g.R1, rc: [g.R6, g.C1], rows: [
          ['同相端分压', 'R3/(R2+R3)', k, 'V/V'], ['同相放大倍数', '1 + R5/R4', A, 'V/V'], ['电压增益', 'R3/(R2+R3)·(1+R5/R4)', k * A, 'V/V'], ['输出零点', 'Vs·R2/(R2+R3)·(1+R5/R4)', g.Vs * g.R2 / (g.R2 + g.R3) * A, 'V']] }; } },
    { id: 'isnH', name: '霍尔电流传感器 (+ 可选分压)', board: '', fields: [
      ['S', '灵敏度', 'mV/A', 40], ['V0', '零电流输出', 'V', 2.5], ['R1', '分压上臂 (0 = 不分压)', 'kΩ', 0], ['R2', '分压下臂', 'kΩ', 10], ['R3', '输出 RC 电阻', 'Ω', 100], ['C1', '输出 RC 电容', 'nF', 1]],
      draw: v => SV.svg(620, 250,
        SV.t(130, 88, '被测电流 i 穿过', 'lv') + SV.box(60, 100, 140, 70, '霍尔传感器', 'S ' + v.S + '，V0 ' + v.V0) + SV.w([200, 135], [240, 135]) +
        SV.R(240, 135, 320, 135, 'R1', v.R1) + SV.dot(340, 135) + SV.w([320, 135], [380, 135]) + SV.R(340, 135, 340, 215, 'R2', v.R2, { left: true }) + SV.gnd(340, 215) +
        SV.R(380, 135, 450, 135, 'R3', v.R3) + SV.w([450, 135], [540, 135]) + SV.dot(490, 135) + SV.Cp(490, 135, 490, 195, 'C1', v.C1) + SV.gnd(490, 195) + SV.port(545, 135, 'ADC')),
      calc: g => { const r = g.R1 > 0 ? g.R2 / (g.R1 + g.R2) : 1;
        return { GVA: g.S * r, Voff: g.V0 * r, rsh: NaN, rc: [g.R3, g.C1], rows: [['分压比', g.R1 > 0 ? 'R2/(R1+R2)' : '无分压 = 1', r, 'V/V']] }; } },
    { id: 'isnP', name: '分流电阻 + MCU 内置 PGA', board: '', fields: [
      ['R1', '采样电阻', 'mΩ', 20], ['R2', '滤波串联电阻 (每边)', 'Ω', 100], ['C1', '差分滤波电容', 'nF', 1], ['Gp', 'PGA 增益', 'V/V', 5], ['Vo', 'PGA 零点电压 (查芯片手册)', 'V', 2.5]],
      draw: v => SV.svg(620, 250,
        SV.t(50, 50, '下桥 N 端') + SV.w([50, 58], [50, 110]) + SV.R(50, 110, 50, 170, 'R1', v.R1, { left: true }) + SV.gnd(50, 170) + SV.dot(50, 110) + SV.dot(50, 170) +
        SV.w([50, 110], [90, 110]) + SV.R(90, 110, 170, 110, 'R2', v.R2) + SV.w([170, 110], [330, 110]) + SV.w([50, 170], [90, 170]) + SV.R(90, 170, 170, 170, 'R2′', v.R2) + SV.w([170, 170], [330, 170]) +
        SV.dot(240, 110) + SV.dot(240, 170) + SV.Cp(240, 110, 240, 170, 'C1', v.C1) + SV.box(330, 90, 180, 100, 'MCU PGA ×' + v.Gp, '零点 ' + v.Vo) + SV.t(342, 114, '+', 'pm', 'start') + SV.t(342, 174, '−', 'pm', 'start') + SV.w([510, 140], [550, 140]) + SV.port(555, 140, 'ADC')),
      calc: g => ({ Gv: g.Gp, Voff: g.Vo, rsh: g.R1, rc: [2 * g.R2, g.C1], rows: [['差分滤波截止', '1/(2π·2R2·C1)', 1 / (TAU * 2 * g.R2 * g.C1) / 1e3, 'kHz']] }) },
  ], common: [['fwGain', '固件 CURRENT_AMP_GAIN', '', 4.545], ['fwOff', '固件 ADC_OFFSET', 'LSB', 2047]] },

  { id: 'vdc', title: '母线电压采样', sub: 'Udc = (Vadc − Vo)·K', topos: [
    { id: 'vdcA', name: '电阻分压 + 电压跟随 (本板)', board: 'U7-B (TLV9062)', fields: [
      ['R1', '上臂电阻 1', 'kΩ', 249, 'R190'], ['R2', '上臂电阻 2', 'kΩ', 249, 'R191'], ['R3', '上臂电阻 3', 'kΩ', 249, 'R192'], ['R4', '上臂电阻 4', 'kΩ', 249, 'R193'],
      ['R5', '下臂电阻', 'kΩ', 10, 'R194'], ['C1', '下臂滤波电容', 'nF', 100, 'C151'], ['R6', '输出 RC 电阻', 'Ω', 100, 'R196'], ['C2', '输出 RC 电容', 'pF', 100, 'C152']],
      draw: v => SV.svg(620, 380,
        SV.t(80, 22, 'HV+ (母线正)') + SV.w([80, 28], [80, 40]) + SV.R(80, 40, 80, 100, 'R1', v.R1) + SV.R(80, 100, 80, 160, 'R2', v.R2) + SV.R(80, 160, 80, 220, 'R3', v.R3) + SV.R(80, 220, 80, 280, 'R4', v.R4) +
        SV.dot(80, 280) + SV.R(80, 280, 80, 340, 'R5', v.R5, { left: true }) + SV.gnd(80, 340) + SV.w([80, 280], [170, 280]) + SV.dot(170, 280) + SV.Cp(170, 280, 170, 340, 'C1', v.C1) + SV.gnd(170, 340) +
        SV.w([170, 280], [260, 280], [260, 270], [300, 270]) + SV.op(300, 250, 'U1', false) + SV.w([300, 230], [280, 230], [280, 195], [380, 195], [380, 250]) +
        SV.w([366, 250], [400, 250]) + SV.dot(380, 250) + SV.R(400, 250, 460, 250, 'R6', v.R6) + SV.w([460, 250], [536, 250]) + SV.dot(500, 250) + SV.Cp(500, 250, 500, 310, 'C2', v.C2) + SV.gnd(500, 310) + SV.port(540, 250, 'ADC')),
      calc: g => { const Rt = g.R1 + g.R2 + g.R3 + g.R4; return { K: (Rt + g.R5) / g.R5, Vo: 0, ng: 1, rth: par(Rt, g.R5), cf: g.C1, rtot: Rt + g.R5, rc: [g.R6, g.C2], rows: [['上臂总电阻', 'R1+R2+R3+R4', Rt / 1e3, 'kΩ'], ['分压比 K', '(R1+…+R4+R5)/R5', (Rt + g.R5) / g.R5, 'V/V']] }; } },
    { id: 'vdcD', name: '电阻分压直接进 ADC', board: '', fields: [['R1', '上臂总电阻', 'kΩ', 1000], ['R2', '下臂电阻', 'kΩ', 10], ['C1', '下臂滤波电容', 'nF', 100]],
      draw: v => SV.svg(620, 260,
        SV.t(80, 22, 'HV+') + SV.w([80, 28], [80, 50]) + SV.R(80, 50, 80, 140, 'R1', v.R1) + SV.dot(80, 140) + SV.R(80, 140, 80, 220, 'R2', v.R2, { left: true }) + SV.gnd(80, 220) +
        SV.w([80, 140], [200, 140]) + SV.dot(200, 140) + SV.Cp(200, 140, 200, 220, 'C1', v.C1) + SV.gnd(200, 220) + SV.w([200, 140], [330, 140]) + SV.port(335, 140, 'ADC')),
      calc: g => ({ K: (g.R1 + g.R2) / g.R2, Vo: 0, rth: par(g.R1, g.R2), cf: g.C1, rtot: g.R1 + g.R2, note: 'ADC 源阻抗 = R1∥R2，需靠 C1 提供采样电荷；C1 越大响应越慢。', rows: [['分压比 K', '(R1+R2)/R2', (g.R1 + g.R2) / g.R2, 'V/V'], ['ADC 源阻抗', 'R1∥R2', par(g.R1, g.R2) / 1e3, 'kΩ']] }) },
    { id: 'vdcN', name: '电阻分压 + 同相放大', board: '', fields: [['R1', '上臂总电阻', 'kΩ', 2000], ['R2', '下臂电阻', 'kΩ', 10], ['C1', '下臂滤波电容', 'nF', 100], ['R3', '增益电阻 (接地)', 'kΩ', 10], ['R4', '反馈电阻', 'kΩ', 10], ['R5', '输出 RC 电阻', 'Ω', 100], ['C2', '输出 RC 电容', 'pF', 100]],
      draw: v => SV.svg(620, 300,
        SV.t(80, 22, 'HV+') + SV.w([80, 28], [80, 50]) + SV.R(80, 50, 80, 150, 'R1', v.R1) + SV.dot(80, 150) + SV.R(80, 150, 80, 230, 'R2', v.R2, { left: true }) + SV.gnd(80, 230) +
        SV.w([80, 150], [170, 150]) + SV.dot(170, 150) + SV.Cp(170, 150, 170, 230, 'C1', v.C1) + SV.gnd(170, 230) + SV.w([170, 150], [300, 150]) +
        SV.op(300, 170, 'U1', true) + SV.w([300, 190], [270, 190]) + SV.dot(270, 190) + SV.R(270, 190, 270, 260, 'R3', v.R3, { left: true }) + SV.gnd(270, 260) +
        SV.w([285, 190], [285, 240]) + SV.dot(285, 190) + SV.R(285, 240, 390, 240, 'R4', v.R4) + SV.w([390, 240], [390, 170]) +
        SV.w([366, 170], [410, 170]) + SV.dot(390, 170) + SV.R(410, 170, 470, 170, 'R5', v.R5) + SV.w([470, 170], [546, 170]) + SV.dot(510, 170) + SV.Cp(510, 170, 510, 230, 'C2', v.C2) + SV.gnd(510, 230) + SV.port(550, 170, 'ADC')),
      calc: g => { const A = 1 + g.R4 / g.R3; return { K: (g.R1 + g.R2) / g.R2 / A, Vo: 0, ng: A, rth: par(g.R1, g.R2), cf: g.C1, rtot: g.R1 + g.R2, rc: [g.R5, g.C2], rows: [['分压比', 'R2/(R1+R2)', g.R2 / (g.R1 + g.R2), 'V/V'], ['放大倍数', '1 + R4/R3', A, 'V/V'], ['总比例 K', '(R1+R2)/R2/(1+R4/R3)', (g.R1 + g.R2) / g.R2 / A, 'V/V']] }; } },
    { id: 'vdcF', name: '差分放大 (HV+ / HV− 两端电阻)', board: '', fields: [['R1', 'HV+ 侧输入电阻', 'kΩ', 2000], ['R2', '同相端接地电阻', 'kΩ', 20], ['R3', 'HV− 侧输入电阻', 'kΩ', 2000], ['R4', '反馈电阻', 'kΩ', 20], ['R5', '输出 RC 电阻', 'Ω', 100], ['C1', '输出 RC 电容', 'pF', 100]],
      draw: v => SV.svg(620, 280,
        SV.t(40, 120, 'HV+', 'lb', 'start') + SV.w([40, 130], [70, 130]) + SV.R(70, 130, 200, 130, 'R1', v.R1, { above: true }) + SV.w([200, 130], [330, 130]) + SV.dot(260, 130) + SV.R(260, 130, 260, 55, 'R2', v.R2) + SV.vs(260, 55, 'AGND') +
        SV.t(40, 160, 'HV−', 'lb', 'start') + SV.w([40, 170], [70, 170]) + SV.R(70, 170, 200, 170, 'R3', v.R3, { below: true }) + SV.w([200, 170], [330, 170]) +
        SV.dot(290, 170) + SV.w([290, 170], [290, 235]) + SV.R(290, 235, 430, 235, 'R4', v.R4) + SV.w([430, 235], [430, 150]) +
        SV.op(330, 150, 'U1', true) + SV.w([396, 150], [450, 150]) + SV.dot(430, 150) + SV.R(450, 150, 510, 150, 'R5', v.R5) + SV.w([510, 150], [556, 150]) + SV.dot(530, 150) + SV.Cp(530, 150, 530, 210, 'C1', v.C1) + SV.gnd(530, 210) + SV.port(560, 150, 'ADC')),
      calc: g => { const o = opLin(g.R1, g.R2, 0, g.R3, g.R4); return { K: 1 / o.a, Vo: 0, ng: 1 + g.R4 / g.R3, rtot: g.R1 + g.R2, rc: [g.R5, g.C1], rows: [
        ['HV+ 系数 a', '(1+R4/R3)·R2/(R1+R2)', o.a, 'V/V'], ['HV− 系数 b', 'R4/R3', o.b, 'V/V'], ['共模增益 a−b', '≈0 表示电阻匹配', o.a - o.b, 'V/V', Math.abs(o.a - o.b) < 0.01 * o.a ? 'ok' : 'warn'], ['比例 K (HV−=0)', '1/a', 1 / o.a, 'V/V']] }; } },
    { id: 'vdcI', name: '电阻分压 + 隔离放大器', board: '', fields: [['R1', '上臂总电阻', 'kΩ', 2000], ['R2', '下臂电阻 (并联隔离放大器输入)', 'kΩ', 10], ['Gi', '隔离放大器增益', 'V/V', 1], ['Vo', '输出零点电压', 'V', 0]],
      draw: v => SV.svg(620, 260,
        SV.t(80, 22, 'HV+') + SV.w([80, 28], [80, 50]) + SV.R(80, 50, 80, 140, 'R1', v.R1) + SV.dot(80, 140) + SV.R(80, 140, 80, 220, 'R2', v.R2, { left: true }) + SV.gnd(80, 220) +
        SV.w([80, 140], [230, 140]) + SV.box(230, 100, 160, 80, '隔离放大器', '×' + v.Gi + '，零点 ' + v.Vo) + SV.w([310, 80], [310, 100]).replace('polyline', 'polyline stroke-dasharray="4 3"') + SV.w([310, 180], [310, 200]).replace('polyline', 'polyline stroke-dasharray="4 3"') + SV.t(310, 215, '隔离', 'lv') + SV.w([390, 140], [480, 140]) + SV.port(485, 140, 'ADC')),
      calc: g => ({ K: (g.R1 + g.R2) / g.R2 / g.Gi, Vo: g.Vo, rtot: g.R1 + g.R2, rows: [['分压比', 'R2/(R1+R2)', g.R2 / (g.R1 + g.R2), 'V/V'], ['总比例 K', '(R1+R2)/R2/Gi', (g.R1 + g.R2) / g.R2 / g.Gi, 'V/V']] }) },
  ], common: [['vov', '软件过压阈值', 'V', 380], ['vuv', '软件欠压阈值', 'V', 160], ['fwVg', '固件 VOLTAGE_GAIN', '', 100.6]] },

  { id: 'oc', title: '硬件过流保护', sub: '动作电流应高于软件过流阈值、并留有余量', topos: [
    { id: 'ocA', name: '母线电阻 + 比较器 + IPM CIN (本板)', board: 'U18-B，PSS15S92F6-AG', fields: [
      ['R1', '母线采样电阻', 'mΩ', 15, 'R88'], ['R2', '基准上拉 (接 Vs)', 'kΩ', 47, 'R200'], ['R3', '基准下拉', 'kΩ', 4.99, 'R199'], ['Vs', '基准电源', 'V', 5, 'AP5V'],
      ['R4', '比较器输入滤波电阻', 'Ω', 470, 'R81'], ['C1', '比较器输入滤波电容', 'pF', 220, 'C81'], ['R5', 'CIN 滤波电阻', 'Ω', 470, 'R82'], ['C2', 'CIN 滤波电容', 'pF', 220, 'C82'], ['Vsc', 'IPM CIN 阈值 Vsc(ref)', 'V', 0.48, 'IPM 手册']],
      draw: v => SV.svg(620, 340,
        SV.t(60, 22, '逆变器 N 母线') + SV.w([60, 28], [60, 120]) + SV.R(60, 120, 60, 190, 'R1', v.R1, { left: true }) + SV.gnd(60, 190) + SV.dot(60, 70) + SV.dot(60, 120) +
        SV.w([60, 70], [110, 70]) + SV.R(110, 70, 200, 70, 'R5', v.R5) + SV.w([200, 70], [320, 70]) + SV.dot(260, 70) + SV.Cp(260, 70, 260, 125, 'C2', v.C2) + SV.gnd(260, 125) + SV.port(325, 70, 'IPM CIN (Vsc ' + v.Vsc + ')') +
        SV.dot(330, 230) + SV.R(330, 230, 330, 160, 'R2', v.R2) + SV.vs(330, 160, 'Vs ' + v.Vs) + SV.R(330, 230, 250, 230, 'R3', v.R3) + SV.gnd(250, 230) + SV.w([330, 230], [360, 230]) +
        SV.w([60, 120], [110, 120], [110, 270], [130, 270]) + SV.R(130, 270, 220, 270, 'R4', v.R4) + SV.w([220, 270], [360, 270]) + SV.dot(280, 270) + SV.Cp(280, 270, 280, 320, 'C1', v.C1) + SV.gnd(280, 320) +
        SV.op(360, 250, 'U1', true) + SV.w([426, 250], [476, 250]) + SV.port(480, 250, 'CIN_FO → POE')),
      calc: g => { const vr = g.Vs * g.R3 / (g.R2 + g.R3); return { trips: [['比较器', vr / g.R1], ['IPM CIN', g.Vsc / g.R1]], rows: [
        ['比较器基准', 'Vs·R3/(R2+R3)', vr, 'V'], ['比较器动作电流', 'Vref/R1', vr / g.R1, 'A'], ['IPM CIN 动作电流', 'Vsc/R1', g.Vsc / g.R1, 'A'],
        ['比较器滤波时间常数', 'R4·C1', g.R4 * g.C1 * 1e9, 'ns'], ['CIN 滤波时间常数', 'R5·C2', g.R5 * g.C2 * 1e9, 'ns']] }; } },
    { id: 'ocC', name: 'IPM CIN 直接检测', board: '', fields: [['R1', '母线采样电阻', 'mΩ', 15], ['R2', 'CIN 滤波电阻', 'Ω', 470], ['C1', 'CIN 滤波电容', 'pF', 220], ['Vsc', 'IPM CIN 阈值', 'V', 0.48]],
      draw: v => SV.svg(620, 220,
        SV.t(60, 22, '逆变器 N 母线') + SV.w([60, 28], [60, 100]) + SV.R(60, 100, 60, 170, 'R1', v.R1, { left: true }) + SV.gnd(60, 170) + SV.dot(60, 100) +
        SV.w([60, 100], [140, 100]) + SV.R(140, 100, 230, 100, 'R2', v.R2) + SV.w([230, 100], [330, 100]) + SV.dot(280, 100) + SV.Cp(280, 100, 280, 160, 'C1', v.C1) + SV.gnd(280, 160) + SV.port(335, 100, 'IPM CIN')),
      calc: g => ({ trips: [['IPM CIN', g.Vsc / g.R1]], rows: [['动作电流', 'Vsc/R1', g.Vsc / g.R1, 'A'], ['滤波时间常数', 'R2·C1', g.R2 * g.C1 * 1e9, 'ns']] }) },
    { id: 'ocW', name: '相电流放大输出 + 窗口比较器', board: '', fields: [['R1', '阈值梯形电阻 上 (接 Vs)', 'kΩ', 2], ['R2', '阈值梯形电阻 中', 'kΩ', 21], ['R3', '阈值梯形电阻 下', 'kΩ', 2], ['Vs', '阈值电源', 'V', 5]],
      draw: v => SV.svg(620, 290,
        SV.vs(90, 30, 'Vs ' + v.Vs) + SV.R(90, 30, 90, 100, 'R1', v.R1, { left: true }) + SV.dot(90, 100) + SV.R(90, 100, 90, 190, 'R2', v.R2, { left: true }) + SV.dot(90, 190) + SV.R(90, 190, 90, 260, 'R3', v.R3, { left: true }) + SV.gnd(90, 260) +
        SV.t(220, 38, '电流放大器输出 Vout', 'lv') + SV.w([220, 48], [220, 230], [300, 230]) + SV.dot(220, 90) + SV.w([220, 90], [300, 90]) +
        SV.w([90, 100], [160, 100], [160, 130], [300, 130]) + SV.w([90, 190], [300, 190]) +
        SV.op(300, 110, 'U1', true) + SV.op(300, 210, 'U2', true) + SV.t(125, 95, 'Vhi', 'lv') + SV.t(125, 185, 'Vlo', 'lv') +
        SV.w([366, 110], [430, 110], [430, 210], [366, 210]) + SV.dot(430, 160) + SV.w([430, 160], [470, 160]) + SV.port(475, 160, '→ 故障 / POE')),
      calc: (g, ctx) => { const S = g.R1 + g.R2 + g.R3, vh = g.Vs * (g.R2 + g.R3) / S, vl = g.Vs * g.R3 / S;
        const i1 = (vh - ctx.Voff) / ctx.GVA, i2 = (vl - ctx.Voff) / ctx.GVA;
        return { trips: [['窗口', Math.min(Math.abs(i1), Math.abs(i2))]], rows: [['上阈值', 'Vs·(R2+R3)/(R1+R2+R3)', vh, 'V'], ['下阈值', 'Vs·R3/(R1+R2+R3)', vl, 'V'],
          ['上阈值对应电流', '(Vhi − Voff)/G', i1, 'A'], ['下阈值对应电流', '(Vlo − Voff)/G', i2, 'A']] }; } },
    { id: 'ocM', name: 'MCU 内置比较器 + DAC 阈值', board: '', fields: [['nb', 'DAC 位数', 'bit', 8], ['code', 'DAC 设定码值', 'LSB', 240], ['Vr', 'DAC 参考电压', 'V', 5]],
      draw: v => SV.svg(620, 200,
        SV.t(150, 40, '电流放大器输出', 'lv') + SV.w([60, 60], [260, 60], [260, 85], [300, 85]) + SV.box(100, 120, 150, 50, 'DAC', v.code + ' / 2^(' + v.nb + ')') + SV.w([250, 145], [270, 145], [270, 125], [300, 125]) +
        SV.op(300, 105, 'CMP', true) + SV.w([366, 105], [420, 105]) + SV.port(425, 105, '→ POE')),
      calc: (g, ctx) => { const vt = g.code / Math.pow(2, g.nb) * g.Vr, i = (vt - ctx.Voff) / ctx.GVA;
        return { trips: [['MCU 比较器', Math.abs(i)]], rows: [['阈值电压', 'code/2^n·Vref', vt, 'V'], ['对应电流', '(Vth − Voff)/G', i, 'A'], ['阈值分辨率', 'Vref/2^n/G', g.Vr / Math.pow(2, g.nb) / Math.abs(ctx.GVA), 'A']] }; } },
  ], common: [['ocMul', '软件过流倍数 (×额定峰值)', '', 1.5]] },

  { id: 'tmp', title: '温度采样', sub: '模块温度 / 散热器温度', topos: [
    { id: 'tmpV', name: 'IPM VOT 线性输出 (本板)', board: 'PSS15S92F6-AG VOT (U9-17)', fields: [['KT', '温度增益', '°C/V', 40, '固件'], ['T0', '温度偏置', '°C', 20, '固件'], ['R1', '输入滤波电阻', 'kΩ', 1, '请按原理图'], ['C1', '输入滤波电容', 'nF', 10, '请按原理图']],
      draw: v => SV.svg(620, 200,
        SV.box(30, 70, 170, 60, 'IPM VOT', 'T = V·' + v.KT + ' − ' + v.T0) + SV.w([200, 100], [200, 100]) + SV.R(200, 100, 290, 100, 'R1', v.R1) + SV.w([290, 100], [400, 100]) + SV.dot(340, 100) + SV.Cp(340, 100, 340, 160, 'C1', v.C1) + SV.gnd(340, 160) + SV.port(405, 100, 'ADC')),
      calc: g => ({ VofT: T => (T + g.T0) / g.KT, lin: true, KT: g.KT, T0: g.T0, rows: [['滤波时间常数', 'R1·C1', g.R1 * g.C1 * 1e3, 'ms']] }) },
    { id: 'tmpL', name: 'NTC 下置 (上拉电阻 + NTC 接地)', board: '', fields: [['Vs', '供电电压', 'V', 5], ['R1', '上拉电阻', 'kΩ', 10], ['R25', 'NTC 25 °C 阻值', 'kΩ', 10], ['B', 'NTC B 值', 'K', 3950], ['C1', '滤波电容', 'nF', 100]],
      draw: v => SV.svg(620, 260,
        SV.vs(100, 30, 'Vs ' + v.Vs) + SV.R(100, 30, 100, 110, 'R1', v.R1) + SV.dot(100, 110) + SV.R(100, 110, 100, 200, 'NTC', v.R25 + ', B' + v.B, { ntc: true }) + SV.gnd(100, 200) +
        SV.w([100, 110], [260, 110]) + SV.dot(260, 110) + SV.Cp(260, 110, 260, 180, 'C1', v.C1) + SV.gnd(260, 180) + SV.w([260, 110], [370, 110]) + SV.port(375, 110, 'ADC')),
      calc: g => { const Rn = T => g.R25 * Math.exp(g.B * (1 / (T + 273.15) - 1 / 298.15));
        return { VofT: T => g.Vs * Rn(T) / (g.R1 + Rn(T)), Rn, tau: par(g.R1, Rn(25)) * g.C1, rows: [['NTC 阻值 @25 °C', 'R25', g.R25 / 1e3, 'kΩ'], ['电压关系', 'V = Vs·Rntc/(R1+Rntc)', NaN, '']] }; } },
    { id: 'tmpH', name: 'NTC 上置 (NTC 接电源 + 下拉电阻)', board: '', fields: [['Vs', '供电电压', 'V', 5], ['R1', '下拉电阻', 'kΩ', 10], ['R25', 'NTC 25 °C 阻值', 'kΩ', 10], ['B', 'NTC B 值', 'K', 3950], ['C1', '滤波电容', 'nF', 100]],
      draw: v => SV.svg(620, 260,
        SV.vs(100, 30, 'Vs ' + v.Vs) + SV.R(100, 30, 100, 110, 'NTC', v.R25 + ', B' + v.B, { ntc: true }) + SV.dot(100, 110) + SV.R(100, 110, 100, 200, 'R1', v.R1) + SV.gnd(100, 200) +
        SV.w([100, 110], [260, 110]) + SV.dot(260, 110) + SV.Cp(260, 110, 260, 180, 'C1', v.C1) + SV.gnd(260, 180) + SV.w([260, 110], [370, 110]) + SV.port(375, 110, 'ADC')),
      calc: g => { const Rn = T => g.R25 * Math.exp(g.B * (1 / (T + 273.15) - 1 / 298.15));
        return { VofT: T => g.Vs * g.R1 / (g.R1 + Rn(T)), Rn, tau: par(g.R1, Rn(25)) * g.C1, rows: [['电压关系', 'V = Vs·R1/(R1+Rntc)', NaN, '']] }; } },
    { id: 'tmpS', name: '线性温度传感器 (MCP9700 / LM35 类)', board: '', fields: [['V0', '0 °C 输出电压', 'V', 0.5], ['Tc', '温度系数', 'mV/°C', 10]],
      draw: v => SV.svg(620, 160,
        SV.box(30, 50, 220, 60, '温度传感器', 'V = ' + v.V0 + ' + ' + v.Tc + '·T') + SV.w([250, 80], [330, 80]) + SV.port(335, 80, 'ADC')),
      calc: g => ({ VofT: T => g.V0 + g.Tc * T, lin: true, KT: 1 / g.Tc, T0: g.V0 / g.Tc, rows: [] }) },
  ], common: [['tot', '过温阈值', '°C', 100]] },

  { id: 'pwm', title: 'PWM 定时器', sub: 'MTU3/MTU4 互补 PWM，三角波载波', common: [
    ['ftim', '定时器时钟 fTIM', '', 80, 'MHz'], ['tdead', '死区时间 td', '', 2, 'µs'], ['tmin', 'IPM 最小脉宽 tmin', '', 0.7, 'µs']] },
];
const HW_IDS = [];
HWG.forEach(gp => { (gp.common || []).forEach(f => HW_IDS.push(f[0])); if (gp.topos) { HW_IDS.push(gp.id + 'Topo'); gp.topos.forEach(tp => tp.fields.forEach(f => HW_IDS.push(tp.id + '_' + f[0]))); } });

function topoVals(tp) {
  const g = {}, v = {};
  tp.fields.forEach(([k, , u]) => { const raw = IN[tp.id + '_' + k]; const x = parseFloat(raw); g[k] = x * (UNIT[u] ?? 1); v[k] = `${raw}${u ? ' ' + u : ''}`; });
  return { g, v };
}
/* ---------- worst-case tolerance analysis ----------
 * Every toleranced component is set to its +tol or −tol limit, over all combinations (2^n corners; above 12 components
 * a fixed pseudo-random set of 4096 corners). For the linear-fractional circuit relations used here the extremes
 * occur at corners, so min/max over the corners is the worst case. */
function tolClass(k, lab, u) {
  if (k === 'Vsc') return 'tolSc';
  if (k === 'R25') return 'tolNtc';
  if (k === 'B' && u === 'K') return 'tolB';
  if (/Ω$/.test(u)) return /采样电阻/.test(lab) ? 'tolSh' : 'tolR';
  if (u === 'V' && /^V[sc]$/.test(k)) return 'tolV';
  return null;
}
function tolCorners(tp, g, tol, fn) {
  const keys = tp.fields.map(f => [f[0], tolClass(f[0], f[1], f[2])]).filter(([k, c]) => c && tol[c] > 0 && isFinite(g[k]) && g[k] !== 0);
  const n = keys.length, N = n <= 12 ? 1 << n : 4096;
  let seed = 12345; const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
  const lo = {}, hi = {};
  for (let i = 0; i < N; i++) {
    const gg = { ...g };
    keys.forEach(([k, c], j) => { const up = n <= 12 ? (i >> j) & 1 : rnd() < 0.5; gg[k] = g[k] * (1 + (up ? 1 : -1) * tol[c] / 100); });
    const r = fn(gg);
    for (const [k, x] of Object.entries(r)) { if (!isFinite(x)) continue; lo[k] = Math.min(lo[k] ?? Infinity, x); hi[k] = Math.max(hi[k] ?? -Infinity, x); }
  }
  return { lo, hi, n, N };
}
const devPct = (lo, hi, nom) => Math.max(Math.abs(lo / nom - 1), Math.abs(hi / nom - 1)) * 100;

function computeHW(P) {
  const gv = id => parseFloat(IN[id]);
  const v = {}, out = {}, extra = {};
  const Inrms = P.In / Math.SQRT2;
  const act = {};
  HWG.forEach(gp => {
    if (!gp.topos) return;
    const tid = IN[gp.id + 'Topo'], tp = gp.topos.find(x => x.id === tid) || gp.topos[0];
    act[gp.id] = { tp, ...topoVals(tp) };
  });
  const tol = {}; ['tolR', 'tolSh', 'tolV', 'tolSc', 'tolNtc', 'tolB'].forEach(k => { tol[k] = Math.max(0, gv(k) || 0); });
  const vos = (gv('opVos') || 0) * 1e-3, vsw = (gv('opSw') || 0) * 1e-3, gbw = (gv('opGbw') || 0) * 1e6;
  const tolInfo = [];
  // ADC
  v.N = Math.round(gv('adcBits')); v.adcMax = Math.pow(2, v.N) - 1; v.vref = gv('adcVref'); v.vLsb = v.vref / v.adcMax;
  out.adc = [['满量程码值', '2^N − 1', v.adcMax, 'LSB'], ['电压分辨率', 'Vref/(2^N−1)', v.vLsb * 1e3, 'mV/LSB'], ['ADC 时钟周期', '1/fAD', 1e3 / gv('adcClk'), 'ns'], ['采样窗口 / 控制周期', 'tAD·fs', gv('tAd') * 1e-6 * P.fs * 100, '%']];
  // current sensing
  { const a = act.isn, r = a.tp.calc(a.g);
    const GVA = r.GVA !== undefined ? r.GVA : r.Gv * r.rsh;
    const Voff = r.Voff, k = Math.abs(GVA);
    v.GVA = GVA; v.Voff = Voff; v.offCode = Voff / v.vLsb; v.aLsb = v.vLsb / k; v.rsh = r.rsh;
    v.gainFw = isFinite(r.rsh) ? Math.abs(r.Gv) : k; v.rshFw = isFinite(r.rsh) ? r.rsh : 1;
    const ilo = Voff / k, ihi = (v.vref - Voff) / k, irng = Math.min(ilo, ihi);
    v.irng = irng;
    const rows = [...r.rows,
      ['电流-电压系数', isFinite(r.rsh) ? 'G·R1' : 'S·分压比', GVA * 1e3, 'mV/A'],
      ['零电流输出', 'Voff', Voff, 'V'],
      ['零电流码值', 'Voff/LSB', v.offCode, 'LSB', Math.abs(v.offCode - gv('fwOff')) <= 2 ? 'ok' : 'warn'],
      ['电流分辨率', 'LSB/|G·R1|', v.aLsb * 1e3, 'mA/LSB'],
      ['可测范围 (输出→0 V)', 'Voff/|G·R1|', GVA > 0 ? -ilo : ilo, 'A'],
      ['可测范围 (输出→Vref)', '(Vref−Voff)/|G·R1|', GVA > 0 ? ihi : -ihi, 'A'],
      ['量程余量', 'min(范围)/Imax', irng / P.Imax, '×', irng / P.Imax >= 1.2 ? 'ok' : irng / P.Imax >= 1 ? 'warn' : 'bad'],
      ['Imax 对应码值变化', 'Imax·|G·R1|/LSB', P.Imax * k / v.vLsb, 'LSB']];
    if (isFinite(r.rsh)) rows.push(['采样电阻压降 @Imax', 'Imax·R1', P.Imax * r.rsh * 1e3, 'mV'], ['采样电阻功耗 @IN', 'IN,rms²·R1 (上限)', Inrms * Inrms * r.rsh, 'W'],
      ['FW 增益误差', '|G|/FW − 1', (Math.abs(r.Gv) / gv('fwGain') - 1) * 100, '%', Math.abs(Math.abs(r.Gv) / gv('fwGain') - 1) < 0.01 ? 'ok' : 'warn']);
    if (r.fb && r.fb[1] > 0) rows.push(['反馈极点', '1/(2π·Rf·Cf)', 1 / (TAU * r.fb[0] * r.fb[1]) / 1e3, 'kHz']);
    if (r.rc && r.rc[0] > 0 && r.rc[1] > 0) { const tau = r.rc[0] * r.rc[1], ts = tau * Math.log(Math.pow(2, v.N + 1));
      rows.push(['输出 RC 截止', '1/(2π·R·C)', 1 / (TAU * tau) / 1e3, 'kHz'], ['RC 建立到 ½LSB', 'τ·ln(2^(N+1))', ts * 1e6, 'µs', ts * 1e6 < gv('tAd') ? 'ok' : 'warn']); }
    // worst case over component tolerances, op-amp offset and output swing
    const tc = tolCorners(a.tp, a.g, tol, gg => { const x = a.tp.calc(gg); const G = x.GVA !== undefined ? x.GVA : x.Gv * x.rsh;
      return { G: Math.abs(G), Voff: x.Voff, rng: Math.min(x.Voff - vsw, v.vref - vsw - x.Voff) / Math.abs(G) }; });
    tolInfo.push(['相电流采样', tc.n, tc.N]);
    const gDev = devPct(tc.lo.G, tc.hi.G, k), ng = r.ng || 1;
    const offV = Math.max(Math.abs(tc.lo.Voff - Voff), Math.abs(tc.hi.Voff - Voff)) + (r.ng ? vos * ng : 0);
    const rngW = tc.lo.rng, ocSwI = gv('ocMul') * P.In;
    v.irngW = rngW; v.gDev = gDev;
    rows.push(
      ['增益最坏误差', `容差角点 ${tc.N} 组`, gDev, '%', gDev < 2 ? 'ok' : gDev < 5 ? 'warn' : 'bad'],
      ['零点最坏偏移', r.ng ? 'ΔVoff(容差) + Vos·噪声增益' : 'ΔVoff(容差)', offV * 1e3, 'mV'],
      ['零点偏移折算电流', 'ΔVoff/|G·R1|', offV / k, 'A', offV / k < 0.05 * P.Imax ? 'ok' : 'warn'],
      ['最坏可测量程', 'min(Voff−Vsw, Vref−Vsw−Voff)/|G|', rngW, 'A', rngW > ocSwI * 1.05 ? 'ok' : rngW > ocSwI ? 'warn' : 'bad']);
    if (r.ng) rows.push(['运放噪声增益', '1 + Rf/Rn', ng, 'V/V'], ['运放闭环带宽', 'GBW / 噪声增益', gbw / ng / 1e3, 'kHz', gbw / ng > 20 * P.fs ? 'ok' : 'warn']);
    out.isn = rows; extra.isn = [r.note, '容差分析：零点偏移可由固件上电标定消除（CURRENT_CFG_OFFSET_CALC_TIME），增益误差不能自动补偿，会直接成为电流和转矩误差。'].filter(Boolean).join(' '); }
  // bus voltage
  { const a = act.vdc, r = a.tp.calc(a.g); v.Kv = r.K; v.Vo = r.Vo || 0; v.vLsbDc = v.vLsb * r.K;
    const code = U => (U / r.K + v.Vo) / v.vLsb;
    const rows = [...r.rows,
      ['满量程母线电压', '(Vref − Vo)·K', (v.vref - v.Vo) * r.K, 'V', (v.vref - v.Vo) * r.K > gv('vov') * 1.05 ? 'ok' : 'bad'],
      ['电压分辨率', 'LSB·K', v.vLsbDc, 'V/LSB'],
      ['额定 Udc 码值', 'Udc/K/LSB', code(P.Vdc), 'LSB'], ['过压阈值码值', 'UOV/K/LSB', code(gv('vov')), 'LSB'], ['欠压阈值码值', 'UUV/K/LSB', code(gv('vuv')), 'LSB'],
      ['FW 误差', 'K/FW − 1', (r.K / gv('fwVg') - 1) * 100, '%', Math.abs(r.K / gv('fwVg') - 1) < 0.005 ? 'ok' : 'warn']];
    if (r.rth && r.cf) rows.push(['分压滤波截止', '1/(2π·(上臂∥下臂)·C)', 1 / (TAU * r.rth * r.cf), 'Hz']);
    if (r.rtot) rows.push(['分压电阻功耗 @UOV', 'UOV²/R总', gv('vov') * gv('vov') / r.rtot * 1e3, 'mW']);
    if (r.rc && r.rc[0] > 0 && r.rc[1] > 0) rows.push(['输出 RC 截止', '1/(2π·R·C)', 1 / (TAU * r.rc[0] * r.rc[1]) / 1e3, 'kHz']);
    const tc = tolCorners(a.tp, a.g, tol, gg => ({ K: a.tp.calc(gg).K }));
    tolInfo.push(['母线电压', tc.n, tc.N]);
    const kDev = devPct(tc.lo.K, tc.hi.K, r.K), vosV = r.ng ? vos * r.ng * r.K : 0, eOV = kDev / 100 * gv('vov') + vosV;
    rows.push(['分压比最坏误差', `容差角点 ${tc.N} 组`, kDev, '%', kDev < 1.5 ? 'ok' : kDev < 3 ? 'warn' : 'bad'],
      ['过压阈值最坏误差', 'ΔK·UOV + Vos·噪声增益·K', eOV, 'V', eOV < 10 ? 'ok' : 'warn']);
    out.vdc = rows; extra.vdc = r.note; }
  // overcurrent
  { const a = act.oc, r = a.tp.calc(a.g, { GVA: v.GVA, Voff: v.Voff }); const ocSw = gv('ocMul') * P.In;
    const rows = [...r.rows, ['软件过流阈值', 'kOC·IN', ocSw, 'A', ocSw < v.irng ? 'ok' : 'bad'], ['软件阈值 / ADC 量程', 'kOC·IN / 范围', ocSw / v.irng * 100, '%', ocSw < v.irng ? 'ok' : 'bad']];
    r.trips.forEach(([n, i]) => rows.push([`${n}动作 / 软件阈值`, 'Itrip / IOC', i / ocSw, '×', i > ocSw * 1.1 ? 'ok' : i > ocSw ? 'warn' : 'bad']));
    const tc = tolCorners(a.tp, a.g, tol, gg => Object.fromEntries(a.tp.calc(gg, { GVA: v.GVA, Voff: v.Voff }).trips.map(([n, i]) => [n, Math.abs(i)])));
    tolInfo.push(['硬件过流', tc.n, tc.N]);
    r.trips.forEach(([n]) => {
      if (!(n in tc.lo)) return;
      rows.push([`${n}动作电流 最小 (容差)`, `容差角点 ${tc.N} 组`, tc.lo[n], 'A', tc.lo[n] > ocSw * 1.05 ? 'ok' : tc.lo[n] > ocSw ? 'warn' : 'bad'],
        [`${n}动作电流 最大 (容差)`, `容差角点 ${tc.N} 组`, tc.hi[n], 'A']);
    });
    if (v.irngW < ocSw) rows.push(['软件阈值超出最坏量程', 'kOC·IN > 最坏可测量程', ocSw - v.irngW, 'A', 'bad']);
    out.oc = rows; }
  // temperature
  { const a = act.tmp, r = a.tp.calc(a.g); const Tot = gv('tot');
    const V = r.VofT, dVdT = (V(Tot + 0.5) - V(Tot - 0.5));
    v.tLsb = Math.abs(v.vLsb / dVdT);
    v.tg = r.lin ? r.KT : 1 / dVdT; v.toff = r.lin ? r.T0 : V(Tot) / dVdT - Tot;
    const rows = [...r.rows.filter(x => isFinite(x[2])),
      ['过温阈值电压', 'V(TOT)', V(Tot), 'V'], ['过温阈值码值', 'V(TOT)/LSB', V(Tot) / v.vLsb, 'LSB'],
      ['TOT 附近温度分辨率', 'LSB / |dV/dT|', v.tLsb, '°C/LSB', v.tLsb < 0.5 ? 'ok' : 'warn'],
      ['线性化 KT (TOT 处)', '1/(dV/dT)', v.tg, '°C/V'], ['线性化 T0 (TOT 处)', 'V·KT − T', v.toff, '°C']];
    if (r.tau) rows.push(['滤波时间常数 (25 °C)', '(R1∥Rntc)·C1', r.tau * 1e3, 'ms']);
    if (!r.lin) {
      // temperature the firmware would read at the true TOT, using the nominal curve to convert back
      const inv = Vx => { let a1 = -55, b1 = 200; const up = V(b1) > V(a1); for (let i = 0; i < 60; i++) { const m = (a1 + b1) / 2; if ((V(m) < Vx) === up) a1 = m; else b1 = m; } return (a1 + b1) / 2; };
      const tc = tolCorners(a.tp, a.g, tol, gg => ({ T: inv(a.tp.calc(gg).VofT(Tot)) }));
      tolInfo.push(['温度', tc.n, tc.N]);
      const dT = Math.max(Math.abs(tc.lo.T - Tot), Math.abs(tc.hi.T - Tot));
      rows.push(['TOT 处最坏温度误差', `容差角点 ${tc.N} 组 (含 R25、B 值)`, dT, '°C', dT < 3 ? 'ok' : dT < 6 ? 'warn' : 'bad']);
    }
    [-20, 0, 25, 50, 75, 100, 125].forEach(T => rows.push([`${T} °C`, r.Rn ? `Rntc = ${fmt(r.Rn(T) / 1e3, 4)} kΩ` : 'V(T)', V(T), 'V']));
    out.tmp = rows; }
  // PWM
  { const cnt = gv('ftim') * 1e6 / (2 * P.fsw), dcnt = gv('tdead') * gv('ftim'); v.cnt = cnt; v.dcnt = dcnt;
    const dVdt = gv('tdead') * 1e-6 * P.fsw * P.Vdc, Dmax = 1 - (gv('tAd') + gv('tdead')) * 1e-6 * P.fsw, Dmin = gv('tmin') * 1e-6 * P.fsw, vEff = (Dmax - Dmin) * P.Vdc / Math.sqrt(3);
    const eN = P.p * P.nN * TAU / 60 * P.psi;
    out.pwm = [['载波周期计数 (半周期)', 'fTIM/(2·fPWM)', cnt, 'cnt'], ['死区计数', 'td·fTIM', dcnt, 'cnt'], ['占空比分辨率', 'log2(计数)', Math.log2(cnt), 'bit'], ['电压分辨率', 'Udc/计数', P.Vdc / cnt, 'V'],
      ['死区平均电压误差', 'td·fPWM·Udc', dVdt, 'V'], ['死区误差 / Rs·IN', 'ΔUdt/(Rs·IN)', dVdt / (P.Rs * P.In), '×', dVdt / (P.Rs * P.In) < 0.5 ? 'ok' : 'warn'],
      ['低边采样限制 Dmax', '1 − (tAD+td)·fPWM', Dmax * 100, '%'], ['最小脉宽 Dmin', 'tmin·fPWM', Dmin * 100, '%'], ['可用相电压幅值 (近似)', '(Dmax−Dmin)·Udc/√3', vEff, 'V'],
      ['线性区上限', 'Udc/√3', P.Vdc / Math.sqrt(3), 'V'], ['额定转速反电势 (相峰值)', 'p·ωm·ψf', eN, 'V', eN < vEff ? 'ok' : 'warn']]; }
  out.tol = tolInfo.map(([n, c, N]) => [`${n}：参与容差的器件`, N === 1 << c ? `全部 ${N} 个角点` : `${N} 个伪随机角点`, c, '个']);
  return { v, out, act, extra };
}
function renesasCfg(R) {
  const { P, cur, spd, obs, hw } = R; const g = id => parseFloat(R.I[id]);
  const f = x => fmt(x, 6), cf = x => cflt(x).replace(/f$/, '');
  const L = P.Lq;
  // equivalent firmware (ωn, ζ) for the chosen design (q axis / speed)
  const wcn = Math.sqrt(cur.q.Ki / L), zc = (cur.q.Kp + P.Rs) / (2 * wcn * L);
  const Ktd = P.Kt / P.kT;
  const wsn = Math.sqrt(spd.Ki * P.kT * Ktd / P.J), zs = spd.Kp * P.kT * Ktd / (2 * wsn * P.J);
  const psiOut = P.tf === 'pow' ? P.psi * Math.sqrt(1.5) : P.psi;
  const L1 = [
    ['/* r_motor_targetmotor_cfg.h */'],
    ['MOTOR_CFG_POLE_PAIRS', String(P.p), 'Number of pole pairs'],
    ['MOTOR_CFG_MAGNETIC_FLUX', cf(psiOut) + 'f', P.tf === 'pow' ? 'ψa = √1.5·ψf [Wb]' : 'ψf [Wb]'],
    ['MOTOR_CFG_RESISTANCE', cf(P.Rs) + 'f', '[ohm]'],
    ['MOTOR_CFG_D_INDUCTANCE', cf(P.Ld) + 'f', '[H]'],
    ['MOTOR_CFG_Q_INDUCTANCE', cf(P.Lq) + 'f', '[H]'],
    ['MOTOR_CFG_ROTOR_INERTIA', cf(P.J) + 'f', '[kgm^2] (电机+负载)'],
    ['MOTOR_CFG_NOMINAL_CURRENT_RMS', cf(P.In / Math.SQRT2) + 'f', '[Arms]'],
    ['MOTOR_CFG_MAX_SPEED_RPM', cf(P.nN) + 'f', '[rpm]'],
    [''], ['/* r_motor_inverter_cfg.h */'],
    ['INVERTER_CFG_SHUNT_RESIST', cf(hw.v.rshFw) + 'f', isFinite(hw.v.rsh) ? '[ohm]' : '霍尔方案：取 1，增益填 V/A'],
    ['INVERTER_CFG_DEADTIME', cf(g('tdead')) + 'f', '[us]'],
    ['INVERTER_CFG_MIN_PULSE_WIDTH', cf(g('tmin')) + 'f', '[us]'],
    ['INVERTER_CFG_VOLTAGE_GAIN', cf(hw.v.Kv) + 'f', '(Rt+Rb)/Rb·Kbuf'],
    ['INVERTER_CFG_TEMP_GAIN', cf(hw.v.tg) + 'f', '线性化于过温阈值处'],
    ['INVERTER_CFG_TEMP_OFFSET', cf(hw.v.toff) + 'f', ''],
    ['INVERTER_CFG_CURRENT_AMP_GAIN', cf(hw.v.gainFw) + 'f', '电压增益 |G|'],
    ['INVERTER_CFG_CURRENT_LIMIT', cf(P.Imax) + 'f', '[A]'],
    ['INVERTER_CFG_OVERVOLTAGE_LIMIT', cf(g('vov')) + 'f', '[V]'],
    ['INVERTER_CFG_UNDERVOLTAGE_LIMIT', cf(g('vuv')) + 'f', '[V]'],
    ['INVERTER_CFG_OVERTEMPERATURE_LIMIT', cf(g('tot')) + 'f', '[degC]'],
    ['INVERTER_CFG_INPUT_V', cf(P.Vdc) + 'f', '[V]'],
    ['INVERTER_CFG_ADC_REF_VOLTAGE', cf(hw.v.vref) + 'f', '[V]'],
    [''], ['/* r_motor_module_cfg.h */'],
    ['MOTOR_MCU_CFG_PWM_TIMER_FREQ', cf(g('ftim')) + 'f', '[MHz]'],
    ['MOTOR_MCU_CFG_CARRIER_FREQ', cf(P.fsw / 1000) + 'f', '[kHz]'],
    ['MOTOR_MCU_CFG_AD_FREQ', cf(g('adcClk')) + 'f', '[MHz]'],
    ['MOTOR_MCU_CFG_AD12BIT_DATA', cf(hw.v.adcMax) + 'f', ''],
    ['MOTOR_MCU_CFG_ADC_OFFSET', '(0x' + Math.floor(hw.v.offCode).toString(16).toUpperCase() + ')', '零电流码值 ' + fmt(hw.v.offCode, 6)],
    ['MOTOR_COMMON_SPEED_LOOP_DIVIDER', '(' + P.Nw + 'U)', ''],
    ['CURRENT_CFG_OMEGA', cf(wcn / TAU) + 'f', cur.m === 'fw' ? '[Hz]' : '[Hz] 由当前设计等效换算'],
    ['CURRENT_CFG_ZETA', cf(zc) + 'f', ''],
    ['CURRENT_CFG_E_OBS_OMEGA', cf(obs.emf.wo / TAU) + 'f', '[Hz]'],
    ['CURRENT_CFG_E_OBS_ZETA', cf(obs.emf.z) + 'f', ''],
    ['CURRENT_CFG_PLL_EST_OMEGA', cf(obs.pll.wn / TAU) + 'f', '[Hz]'],
    ['CURRENT_CFG_PLL_EST_ZETA', cf(obs.pll.z) + 'f', ''],
    ['SPEED_CFG_OMEGA', cf(wsn / TAU) + 'f', spd.m === 'fw' ? '[Hz]' : '[Hz] 由当前设计等效换算'],
    ['SPEED_CFG_ZETA', cf(zs) + 'f', ''],
    ['SPEED_CFG_LPF_OMEGA', cf(P.ff) + 'f', '[Hz]'],
  ];
  let s = `/* 瑞萨 RX24T 电机库配置宏 — 由 PMSM 控制参数整定工具生成 (${new Date().toISOString().slice(0, 10)})\n * 按需替换对应头文件中的同名宏；修改前请核对。\n * 固件在线计算: Kp = 2ζω·L − R, ki = ω²·L·Ts (电流环); Kp = 2ζω·J/(p·ψa), ki = ω²·J·Tω/(p·ψa) (速度环)\n */\n`;
  L1.forEach(r => { if (r.length === 1) s += r[0] + '\n'; else s += `#define     ${r[0].padEnd(36)}(${r[1].replace(/^\(|\)$/g, '')})${r[2] ? '   /* ' + r[2] + ' */' : ''}\n`; });
  s += `\n/* 固件实际 PI 增益 (运行时由上述 ω、ζ 计算):\n * id: kp = ${f(2 * zc * wcn * P.Ld - P.Rs)}, ki = ${f(wcn * wcn * P.Ld * P.Ts)}\n * iq: kp = ${f(2 * zc * wcn * P.Lq - P.Rs)}, ki = ${f(wcn * wcn * P.Lq * P.Ts)}\n * speed: kp = ${f(2 * zs * wsn * P.J / Ktd)}, ki = ${f(wsn * wsn * P.J / Ktd * P.Tw)}\n * i = (ADC_OFFSET − ad) × ${f(hw.v.aLsb)} A/LSB;  Udc = ad × ${f(hw.v.vLsbDc)} V/LSB\n */\n`;
  return s;
}


/* ============ DEFAULTS & ENTRY POINT ============ */
const BASE_DEFAULTS = {"motorType": "ipm", "tf": "pow", "zcc": "1", "zsc": "1", "Rs": "1.34", "Ld": "8.2", "Lq": "11.23", "pp": "3", "fluxMode": "psia", "fluxVal": "0.086", "Jm": "3", "JL": "0", "Bv": "0.0001", "Irated": "14.14", "Imax": "20", "nRated": "3000", "nMax": "3500", "Vdc": "310", "fsw": "5000", "upd": "single", "nd": "1", "disc": "be", "Nw": "5", "ff": "100", "Np": "10", "enc": "131072", "Vb": "", "Ib": "", "nb": "", "fcc": "20", "curMethod": "fw", "curStep": "4", "spdMethod": "fw", "fsc": "2", "Nz": "4", "aSO": "2.4", "ktMode": "id0", "spdStep": "100", "loadPct": "50", "fpc": "0.5", "vff": "0", "posStep": "30", "fpll": "25", "zpll": "1", "femf": "250", "zemf": "1", "smoEta": "1.5", "smoLpf": "1000", "flo": "50", "zlo": "0.707", "piImpl": "rn", "ilimC": "", "ilimS": "", "spdSrc": "pll", "tgtPmC": "60", "tgtPmS": "50", "tgtGm": "6", "rbL": "30", "rbRlo": "0.9", "rbRhi": "1.5", "rbJlo": "0.5", "rbJhi": "3", "rbKt": "10", "tuneMode": "worst"};
const DEFAULTS = (() => {
  const d = { ...BASE_DEFAULTS };
  HWG.forEach(gp => {
    (gp.common || []).forEach(f => { d[f[0]] = String(f[3]); });
    if (gp.topos) { d[gp.id + 'Topo'] = gp.topos[0].id; gp.topos.forEach(tp => tp.fields.forEach(f => { d[tp.id + '_' + f[0]] = String(f[3]); })); }
  });
  return d;
})();
const FIELD_IDS = Object.keys(DEFAULTS);

/* Runs every design and analysis for one set of inputs. Inputs missing from `inputs` take their defaults. */
function computeAll(inputs) {
  IN = { ...DEFAULTS, ...inputs };
  const I = IN;
  const P = readP();
  const bad = validP(P);
  if (bad.length) return { I, P, bad };
  const cur = designCurrent(P);
  const spd = designSpeed(P, cur);
  const pos = designPos(P, cur, spd);
  const obs = designObs(P);
  const R = { I, P, cur, spd, pos, obs };
  R.hw = computeHW(P);
  return R;
}

/* ============ CODE EXPORT ============ */
function exportValues(R) {
  const { P, cur, spd, pos, obs, hw } = R, k = P.kT;
  return [
    ['基本参数', [
      ['MOTOR_RS', P.Rs, 'Ω'], ['MOTOR_LD', P.Ld, 'H'], ['MOTOR_LQ', P.Lq, 'H'], ['MOTOR_PSI_F', P.psi, 'Wb'], ['MOTOR_POLE_PAIRS', P.p, ''],
      ['MOTOR_J', P.J, 'kg·m²'], ['MOTOR_B', P.B, 'N·m·s/rad'], ['MOTOR_KT', P.Kt, 'N·m/A'], ['VDC', P.Vdc, 'V'], ['V_MAX', P.Vmax * k, 'V, dq 电压上限'], ['I_MAX', P.Imax * k, 'A, dq 电流上限'],
      ['CTRL_TS', P.Ts, 's, 电流环周期'], ['CTRL_DELAY', P.Nd, '拍, 采样到占空比生效'], ['SPD_TS', P.Tw, 's, 速度环周期'], ['POS_TS', P.Tp, 's, 位置环周期']]],
    [`电流环 (fcc=${fmt(cur.fcc, 4)} Hz, ${{ fw: '瑞萨极点配置', zd: 'z域对消', imc: 'IMC' }[cur.m]}, ${DISC_NAME[P.disc]})`, [
      ['ID_KP', cur.d.Kp, 'V/A'], ['ID_KI', cur.d.Ki, 'V/(A·s)'], ['ID_KI_TS', cur.d.Ki * P.Ts, 'Ki·Ts'], ['ID_B0', cur.d.b0, '增量式'], ['ID_B1', cur.d.b1, '增量式'],
      ['IQ_KP', cur.q.Kp, 'V/A'], ['IQ_KI', cur.q.Ki, 'V/(A·s)'], ['IQ_KI_TS', cur.q.Ki * P.Ts, 'Ki·Ts'], ['IQ_B0', cur.q.b0, '增量式'], ['IQ_B1', cur.q.b1, '增量式'],
      ['ID_KP_PU', cur.d.Kp * P.Ib / P.Vb, '标幺'], ['ID_KI_TS_PU', cur.d.Ki * P.Ts * P.Ib / P.Vb, '标幺'], ['IQ_KP_PU', cur.q.Kp * P.Ib / P.Vb, '标幺'], ['IQ_KI_TS_PU', cur.q.Ki * P.Ts * P.Ib / P.Vb, '标幺']]],
    [`速度环 (fsc=${fmt(spd.fsc, 4)} Hz, PM=${fmt(spd.sw.pm, 3)}°)`, [
      ['SPD_KP', spd.Kp * k, 'A/(rad/s) 机械'], ['SPD_KI', spd.Ki * k, 'A/rad'], ['SPD_KI_TS', spd.Ki * k * P.Tw, 'Ki·Tω'], ['SPD_B0', spd.b0 * k, '增量式'], ['SPD_B1', spd.b1 * k, '增量式'],
      ['SPD_FILT_ALPHA', spd.alpha, '一阶IIR系数'], ['SPD_KP_PU', spd.Kp * P.wmb / P.Ib, '标幺'], ['SPD_KI_TS_PU', spd.Ki * P.Tw * P.wmb / P.Ib, '标幺']]],
    [`位置环 (fpc=${fmt(pos.fpc, 4)} Hz)`, [
      ['POS_KP', pos.Kpp, '1/s'], ['POS_KVFF', pos.vff, ''], ['POS_KAFF', pos.Kaff, 'A/(rad/s²)'], ['POS_KP_PU', pos.Kpp * TAU / P.wmb, '标幺']]],
    ['观测器', [
      ['PLL_KP', obs.pll.Kp, '1/s'], ['PLL_KI_TS', obs.pll.Ki * P.Ts, ''],
      ['EMF_KP', obs.emf.Kp, 'Ω'], ['EMF_KI_TS', obs.emf.Ki * P.Ts, 'Ω'], ['EMF_TS_LQ', P.Ts / P.Lq, 'Ts/Lq'],
      ['SMO_K', obs.smo.k * k, 'V'], ['SMO_LPF_ALPHA', obs.smo.alpha, ''],
      ['TLOBS_L1_TS', obs.lo.l1 * P.Tw, ''], ['TLOBS_L2_TS', obs.lo.l2 * P.Tw, ''], ['TLOBS_KCOMP', obs.lo.Kcomp * k, 'A/(N·m)']]],
    ['硬件采样换算', [
      ['ADC_FULL_SCALE', hw.v.adcMax, 'LSB'], ['ADC_VREF', hw.v.vref, 'V'], ['ADC_I_OFFSET', Math.floor(hw.v.offCode), 'LSB, 零电流码值'],
      ['ADC_A_PER_LSB', hw.v.aLsb, 'A/LSB, i = (OFFSET − ad)·k'], ['ADC_VDC_PER_LSB', hw.v.vLsbDc, 'V/LSB'], ['ADC_TEMP_PER_LSB', hw.v.tLsb, '°C/LSB'],
      ['PWM_PERIOD_CNT', Math.round(hw.v.cnt), '三角波半周期计数'], ['PWM_DEADTIME_CNT', Math.round(hw.v.dcnt), '']]],
  ];
}
function exportCode(R, kind, today) {
  if (kind === 'rn') return renesasCfg(R);
  const { P, cur } = R, f = x => fmt(x, 7), vals = exportValues(R), now = today || new Date().toISOString().slice(0, 10);
  if (kind === 'c') {
    let s = `/* PMSM 控制参数 — 由 PMSM 控制参数整定工具生成 (${now})\n * 单位: SI, 速度为机械 rad/s, dq 量按${P.tf === 'pow' ? '等功率变换 (|idq| = √3·Irms)' : '等幅值变换 (|idq| = Ipk)'}\n * 增量式 PI: u[k] = sat(u[k-1] + B0*e[k] + B1*e[k-1])\n */\n#ifndef PMSM_PARAMS_H\n#define PMSM_PARAMS_H\n\n`;
    vals.forEach(([h, rows]) => { s += `/* ---- ${h} ---- */\n`; rows.forEach(([n, v, u]) => { s += `#define ${n.padEnd(18)} (${cflt(v)})${u ? `   /* ${u} */` : ''}\n`; }); s += '\n'; });
    return s + `#endif /* PMSM_PARAMS_H */\n`;
  }
  if (kind === 'm') {
    const d = cur.d, m = Math.floor(P.Nd + 1e-9);
    let s = `%% PMSM 控制参数 — 由 PMSM 控制参数整定工具生成 (${now})\n% 可直接用于 Simulink: Discrete PID 块 (Parallel, ${DISC_NAME[P.disc]}) 填 Kp / Ki, Sample time 填对应 Ts\n\n`;
    vals.forEach(([h, rows]) => { s += `%% ${h}\n`; rows.forEach(([n, v, u]) => { s += `${n} = ${f(v)};${u ? ` % ${u}` : ''}\n`; }); s += '\n'; });
    s += `%% d 轴离散开环 (需 Control System Toolbox)\n% 对象: ZOH 等效 + 延时 Nd = ${fmt(P.Nd, 4)} 拍 (修正 z 变换)\nz = tf('z', CTRL_TS);\nad = ${f(d.a)}; g0 = ${f(d.g0)}; g1 = ${f(d.g1)};\nGd = z^-${m}*(g0*z^-1 + g1*z^-2)/(1 - ad*z^-1);\nCd = (ID_B0 + ID_B1*z^-1)/(1 - z^-1);\nLd_open = Cd*Gd;  % margin(Ld_open); step(feedback(Ld_open,1));\n`;
    return s;
  }
  const inputs = {}; FIELD_IDS.forEach(id => { inputs[id] = R.I[id]; });
  const res = {}; vals.forEach(([, rows]) => rows.forEach(([n, v]) => { res[n] = +f(v); }));
  return JSON.stringify({ inputs, results: res }, null, 2);
}

/* ============ FIRMWARE CONFIG IMPORT ============ */
/* Minimal C preprocessor for configuration headers: strips comments, follows #if/#ifdef/#elif/#else/#endif,
 * collects object-like #define values and evaluates them as arithmetic expressions. No eval() is used. */
const CFG_PREDEF = { MTR_SQRT_2: Math.SQRT2, MTR_SQRT_3: Math.sqrt(3), MTR_TWOPI: 2 * Math.PI, MTR_PI: Math.PI, MTR_ENABLE: 1, MTR_DISABLE: 0 };

// recursive-descent evaluator: || && ! == != < > <= >= + - * / and parentheses
function cArith(src) {
  const s = src; let i = 0;
  const ws = () => { while (i < s.length && /\s/.test(s[i])) i++; };
  const peek = t => { ws(); return s.startsWith(t, i); };
  const eat = t => { if (peek(t)) { i += t.length; return true; } return false; };
  function num() { ws(); const m = /^(\d*\.?\d+(?:[eE][+-]?\d+)?)/.exec(s.slice(i)); if (!m) throw new Error('num'); i += m[1].length; return parseFloat(m[1]); }
  function un() { if (eat('!')) return un() ? 0 : 1; if (eat('-')) return -un(); if (eat('+')) return un(); if (eat('(')) { const v = or(); if (!eat(')')) throw new Error(')'); return v; } return num(); }
  function mul() { let v = un(); for (;;) { if (eat('*')) v *= un(); else if (eat('/')) v /= un(); else return v; } }
  function add() { let v = mul(); for (;;) { if (eat('+')) v += mul(); else if (eat('-')) v -= mul(); else return v; } }
  function rel() { let v = add(); for (;;) { if (eat('<=')) v = +(v <= add()); else if (eat('>=')) v = +(v >= add()); else if (eat('<')) v = +(v < add()); else if (eat('>')) v = +(v > add()); else return v; } }
  function eq() { let v = rel(); for (;;) { if (eat('==')) v = +(v === rel()); else if (eat('!=')) v = +(v !== rel()); else return v; } }
  function and() { let v = eq(); while (eat('&&')) { const r = eq(); v = +(!!v && !!r); } return v; }
  function or() { let v = and(); while (eat('||')) { const r = and(); v = +(!!v || !!r); } return v; }
  try { const v = or(); ws(); return i === s.length ? v : NaN; } catch (e) { return NaN; }
}
function cfgNormalize(expr) {
  return expr
    .replace(/\b0[xX][0-9a-fA-F]+[uUlL]*\b/g, h => String(parseInt(h.replace(/[uUlL]+$/, ''), 16)))
    .replace(/(?<![\w.])(\d*\.?\d+(?:[eE][+-]?\d+)?)[fFuUlL]+\b/g, '$1')
    .replace(/\(\s*(?:u?int\d+_t|float|double|int|unsigned|long|short|char)\s*\)/g, '');
}
function cfgEval(expr, defs, vals, depth, forIf) {
  if (depth > 32) return NaN;
  let s = cfgNormalize(expr);
  if (forIf) s = s.replace(/defined\s*\(\s*(\w+)\s*\)|defined\s+(\w+)/g, (m, a, b) => ((a || b) in defs || (a || b) in CFG_PREDEF ? '1' : '0'));
  let missing = false;
  s = s.replace(/(?<![\w.])[A-Za-z_]\w*/g, id => {   // not inside numbers such as 1.2e-4
    // macros defined in the headers take precedence over the built-in library constants
    if (id in vals) return '(' + vals[id] + ')';
    if (id in defs) { const v = cfgEval(defs[id], defs, vals, depth + 1, false); if (isFinite(v)) { vals[id] = v; return '(' + v + ')'; } }
    if (id in CFG_PREDEF) return '(' + CFG_PREDEF[id] + ')';
    if (forIf) return '0';
    missing = true; return '0';
  });
  if (missing || !s.trim()) return NaN;
  return cArith(s);
}
function parseCfg(text) {
  const s = text.replace(/\/\*[\s\S]*?\*\//g, ' ').replace(/\/\/[^\n]*/g, '').replace(/\\\r?\n/g, ' ');
  const defs = {}, stack = [];
  const active = () => stack.every(x => x.active);
  for (const line of s.split(/\r?\n/)) {
    const m = /^\s*#\s*(\w+)\s*(.*)$/.exec(line); if (!m) continue;
    const d = m[1], rest = m[2].trim();
    if (d === 'ifdef' || d === 'ifndef') { const has = rest.split(/\s+/)[0] in defs; const c = d === 'ifdef' ? has : !has; stack.push({ active: c, taken: c }); }
    else if (d === 'if') { const c = active() && !!cfgEval(rest, defs, {}, 0, true); stack.push({ active: c, taken: c }); }
    else if (d === 'elif') { const t = stack[stack.length - 1]; if (!t) continue; if (t.taken) t.active = false; else { const c = !!cfgEval(rest, defs, {}, 0, true); t.active = c; t.taken = c; } }
    else if (d === 'else') { const t = stack[stack.length - 1]; if (t) { t.active = !t.taken; t.taken = true; } }
    else if (d === 'endif') stack.pop();
    else if (d === 'define' && active()) { const mm = /^(\w+)(\([^)]*\))?\s*(.*)$/.exec(rest); if (mm && !mm[2]) defs[mm[1]] = mm[3].trim(); }
    else if (d === 'undef' && active()) delete defs[rest.split(/\s+/)[0]];
  }
  const vals = {}, unresolved = [];
  Object.keys(defs).forEach(k => { if (!defs[k]) return; const v = cfgEval(defs[k], defs, vals, 0, false); if (isFinite(v)) vals[k] = v; else unresolved.push(k); });
  return { defs, vals, unresolved };
}
// order in which headers are concatenated: later files may use macros from earlier ones
function cfgSortFiles(files) {
  const rank = n => { const s = n.toLowerCase(); return s.includes('common') ? 0 : s.includes('targetmotor') ? 1 : s.includes('inverter') ? 2 : s.includes('module') ? 3 : 4; };
  return files.slice().sort((a, b) => rank(a.name) - rank(b.name));
}

/* macro → input field. conv converts the firmware value to the field unit; set lists extra fields that are implied. */
const CFG_MAP = [
  ['MOTOR_CFG_POLE_PAIRS', 'pp'],
  ['MOTOR_CFG_MAGNETIC_FLUX', 'fluxVal', null, { fluxMode: 'psia', tf: 'pow' }, 'ψa，按等功率变换'],
  ['MOTOR_CFG_RESISTANCE', 'Rs'],
  ['MOTOR_CFG_D_INDUCTANCE', 'Ld', v => v * 1e3, null, 'H → mH'],
  ['MOTOR_CFG_Q_INDUCTANCE', 'Lq', v => v * 1e3, null, 'H → mH'],
  ['MOTOR_CFG_ROTOR_INERTIA', 'Jm', v => v * 1e4, { JL: '0' }, 'kg·m² → kg·cm²，负载惯量置 0'],
  ['MOTOR_CFG_NOMINAL_CURRENT_RMS', 'Irated', v => v * Math.SQRT2, null, 'Arms → 峰值'],
  ['MOTOR_CFG_MAX_SPEED_RPM', 'nRated'],
  ['SPEED_CFG_SPEED_LIMIT_RPM', 'nMax'],
  ['INVERTER_CFG_CURRENT_LIMIT', 'Imax'],
  ['INVERTER_CFG_INPUT_V', 'Vdc'],
  ['MOTOR_MCU_CFG_CARRIER_FREQ', 'fsw', v => v * 1000, { upd: 'single' }, 'kHz → Hz'],
  ['MOTOR_COMMON_SPEED_LOOP_DIVIDER', 'Nw'],
  ['CURRENT_CFG_OMEGA', 'fcc', null, { curMethod: 'fw', piImpl: 'rn', disc: 'be' }],
  ['CURRENT_CFG_ZETA', 'zcc'],
  ['CURRENT_CFG_PI_INTEGRAL_LIMIT_VQ', 'ilimC'],
  ['SPEED_CFG_OMEGA', 'fsc', null, { spdMethod: 'fw' }],
  ['SPEED_CFG_ZETA', 'zsc'],
  ['SPEED_CFG_LPF_OMEGA', 'ff', null, null, 'Hz'],
  ['MOTOR_COMMON_PI_INTEGRAL_LIMIT_IQ', 'ilimS'],
  ['CURRENT_CFG_E_OBS_OMEGA', 'femf'], ['CURRENT_CFG_E_OBS_ZETA', 'zemf'],
  ['CURRENT_CFG_PLL_EST_OMEGA', 'fpll'], ['CURRENT_CFG_PLL_EST_ZETA', 'zpll'],
  ['MOTOR_COMMON_CFG_OVERCURRENT_MARGIN_MULT', 'ocMul'],
  ['INVERTER_CFG_CURRENT_AMP_GAIN', 'fwGain', null, null, '用于与原理图增益对照'],
  ['MOTOR_MCU_CFG_ADC_OFFSET', 'fwOff'],
  ['INVERTER_CFG_VOLTAGE_GAIN', 'fwVg'],
  ['INVERTER_CFG_ADC_REF_VOLTAGE', 'adcVref'],
  ['MOTOR_MCU_CFG_AD12BIT_DATA', 'adcBits', v => Math.round(Math.log2(v + 1)), null, '满量程码值 → 位数'],
  ['MOTOR_MCU_CFG_AD_FREQ', 'adcClk'],
  ['MOTOR_MCU_CFG_PWM_TIMER_FREQ', 'ftim'],
  ['INVERTER_CFG_DEADTIME', 'tdead'],
  ['INVERTER_CFG_MIN_PULSE_WIDTH', 'tmin'],
  ['INVERTER_CFG_OVERVOLTAGE_LIMIT', 'vov'],
  ['INVERTER_CFG_UNDERVOLTAGE_LIMIT', 'vuv'],
  ['INVERTER_CFG_OVERTEMPERATURE_LIMIT', 'tot'],
  ['INVERTER_CFG_TEMP_GAIN', 'tmpV_KT', null, { tmpTopo: 'tmpV' }],
  ['INVERTER_CFG_TEMP_OFFSET', 'tmpV_T0'],
];
const num2str = v => String(+(+v).toPrecision(7));
/* Applies parsed macro values onto `base` inputs. Returns the new inputs, what was applied and what needs attention. */
function cfgToInputs(vals, base) {
  const inputs = { ...base }, applied = [], notes = [];
  CFG_MAP.forEach(([mac, field, conv, set, note]) => {
    if (!(mac in vals)) return;
    const v = conv ? conv(vals[mac]) : vals[mac];
    inputs[field] = num2str(v);
    if (set) Object.assign(inputs, set);
    applied.push([mac, field, vals[mac], inputs[field], note || '']);
  });
  // shunt resistor goes to R1 of the active current-sensing topology when it is a shunt circuit
  if ('INVERTER_CFG_SHUNT_RESIST' in vals) {
    const tp = HWG.find(g => g.id === 'isn').topos.find(t => t.id === inputs.isnTopo);
    const f = tp && tp.fields.find(x => x[0] === 'R1' && x[2] === 'mΩ');
    if (f) { inputs[tp.id + '_R1'] = num2str(vals.INVERTER_CFG_SHUNT_RESIST * 1e3); applied.push(['INVERTER_CFG_SHUNT_RESIST', tp.id + '_R1', vals.INVERTER_CFG_SHUNT_RESIST, inputs[tp.id + '_R1'], 'Ω → mΩ']); }
    else notes.push('当前相电流采样形式没有分流电阻，INVERTER_CFG_SHUNT_RESIST 未导入。');
  }
  if (vals.MOTOR_MCU_CFG_INTR_DECIMATION > 0) notes.push(`MOTOR_MCU_CFG_INTR_DECIMATION = ${vals.MOTOR_MCU_CFG_INTR_DECIMATION}：控制频率为载波的 1/${vals.MOTOR_MCU_CFG_INTR_DECIMATION + 1}，请按实际控制周期核对 fpwm。`);
  if (!('MOTOR_CFG_MAGNETIC_FLUX' in vals)) notes.push('未找到 MOTOR_CFG_MAGNETIC_FLUX：坐标变换与磁链保持原值。');
  if (vals.MOTOR_COMMON_CFG_LOOP_MODE !== undefined && vals.MOTOR_LOOP_SPEED !== undefined && vals.MOTOR_COMMON_CFG_LOOP_MODE !== vals.MOTOR_LOOP_SPEED) notes.push('固件不是速度环模式，速度环参数仅供参考。');
  return { inputs, applied, notes };
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    DEFAULTS, FIELD_IDS, HWG, computeAll, exportCode, exportValues, renesasCfg, parseCfg, cfgToInputs, cfgSortFiles, cArith,
    piCoef, plantZ, makePI, roots, sweep, fmt, cflt, TAU, autoTune, tradeoff, robustness, currentLoop, speedLoop, pllResp,
  };
}
