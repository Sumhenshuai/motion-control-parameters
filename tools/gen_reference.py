"""Independent reference values for the current-loop analysis (python-control).

Reads tests/fixtures/cases.json (written by tools/make_cases.js) and rebuilds each loop from its physical
parameters only: ZOH-discretised 1/(Ls+R), an integer sample delay z^-Nd and the PI in the chosen discretisation.
Writes tests/fixtures/reference.json. Run:  node tools/make_cases.js && python3 tools/gen_reference.py
"""
import json, math, pathlib
import numpy as np
import control as ct

root = pathlib.Path(__file__).resolve().parent.parent / 'tests' / 'fixtures'
cases = json.loads((root / 'cases.json').read_text())
out = []
for c in cases:
    Ts, nd = c['Ts'], int(round(c['Nd']))
    assert abs(c['Nd'] - nd) < 1e-12, 'python-control reference only covers integer delays'
    plant = ct.c2d(ct.tf([1.0], [c['L'], c['Rs']]), Ts, method='zoh')
    delay = ct.tf([1.0], [1.0] + [0.0] * nd, Ts)
    z = ct.tf([1.0, 0.0], [1.0], Ts)
    Kp, Ki = c['Kp'], c['Ki']
    if c['disc'] == 'fe':
        pi = Kp + Ki * Ts / (z - 1)
    elif c['disc'] == 'tu':
        pi = Kp + Ki * Ts / 2 * (z + 1) / (z - 1)
    else:
        pi = Kp + Ki * Ts * z / (z - 1)
    if c['name'].startswith('z-domain'):          # z-domain design is defined by its coefficients, not Kp/Ki
        pi = ct.tf([c['b0'], c['b1']], [1.0, -1.0], Ts)
    L_full = pi * plant * delay                    # no cancellation: keeps near pole-zero pairs (slow disturbance modes)
    L = ct.minreal(L_full, verbose=False)
    gm, pm, wcg, wcp = ct.margin(L)
    cl = ct.feedback(L_full, 1)
    poles = ct.poles(cl)
    out.append({**c, 'pm': float(pm), 'gm_db': float(20 * math.log10(gm)) if np.isfinite(gm) and gm > 1 else None,
                'fc_hz': float(wcp / (2 * math.pi)), 'max_pole': float(max(abs(p) for p in poles))})
(root / 'reference.json').write_text(json.dumps(out, indent=1))
print(len(out), 'reference cases written')
for o in out:
    print(f"{o['name'][:60]:60s} PM {o['pm']:7.3f}  GM {o['gm_db'] if o['gm_db'] is None else round(o['gm_db'],3)}  fc {o['fc_hz']:8.3f}  |p| {o['max_pole']:.5f}")
