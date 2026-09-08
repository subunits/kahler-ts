"""
kerr_isco.py — Kerr ISCO geodesic integrator
Spin parameter a/M = 0.7
Units: G = c = M = 1

Integrates the Kerr geodesic equations in Boyer-Lindquist coordinates
for a prograde circular orbit at the ISCO radius, with small oscillatory
perturbations added so all four channels carry non-trivial signal.

Output: CSV files written to the current directory.
  isco_t_coord.csv
  isco_r_coord.csv
  isco_theta_coord.csv
  isco_phi_coord.csv
  isco_all.csv

Usage:
  python3 kerr_isco.py
"""

import math
import csv
import os

# ── Kerr parameters ──────────────────────────────────────────────────
M = 1.0
a = 0.7  # spin a/M = 0.7

# ISCO radius for prograde orbit (exact formula, Bardeen 1972)
Z1 = 1 + (1 - a**2)**(1/3) * ((1+a)**(1/3) + (1-a)**(1/3))
Z2 = math.sqrt(3*a**2 + Z1**2)
r_isco = M * (3 + Z2 - math.sqrt((3 - Z1)*(3 + Z1 + 2*Z2)))
print(f"ISCO radius for a/M={a}: r_ISCO = {r_isco:.6f} M")

# Specific energy and angular momentum for circular geodesic
r = r_isco
E = (r**2 - 2*M*r + a*math.sqrt(M*r)) / \
    (r * math.sqrt(r**2 - 3*M*r + 2*a*math.sqrt(M*r)))
L = math.sqrt(M*r) * (r**2 - 2*a*math.sqrt(M*r) + a**2) / \
    (r * math.sqrt(r**2 - 3*M*r + 2*a*math.sqrt(M*r)))
print(f"Energy E = {E:.6f}")
print(f"Ang. mom. L = {L:.6f}")

# ── Geodesic equations (RHS) ─────────────────────────────────────────
def kerr_rhs(lam, state, E, L, M, a):
    t, r, theta, phi = state
    sin_t = math.sin(theta)
    Sigma  = r**2 + a**2 * math.cos(theta)**2
    Delta  = r**2 - 2*M*r + a**2
    dt_dlam  = (E*(r**2+a**2) - a*L) / (Sigma*Delta) * (r**2+a**2) \
               - a*(a*E*sin_t**2 - L) / Sigma
    dr_dlam  = 0.0
    dth_dlam = 0.0
    dph_dlam = (L / (Sigma*sin_t**2) - a*E/Sigma) \
               + a*(E*(r**2+a**2) - a*L) / (Sigma*Delta)
    return [dt_dlam, dr_dlam, dth_dlam, dph_dlam]

# ── RK4 integrator ───────────────────────────────────────────────────
def rk4_step(f, lam, state, h, *args):
    k1 = f(lam,       state,                                *args)
    k2 = f(lam+h/2,   [s+h/2*k for s,k in zip(state,k1)],  *args)
    k3 = f(lam+h/2,   [s+h/2*k for s,k in zip(state,k2)],  *args)
    k4 = f(lam+h,     [s+h*k   for s,k in zip(state,k3)],  *args)
    return [s+h/6*(a+2*b+2*c+d)
            for s,a,b,c,d in zip(state,k1,k2,k3,k4)]

# ── Integrate ────────────────────────────────────────────────────────
state = [0.0, r_isco, math.pi/2, 0.0]
h = 0.5
N = 300

rows = []
for i in range(N):
    lam   = i * h
    eps_r  = 0.02 * r_isco * math.sin(2*math.pi*i/40)
    eps_th = 0.015 * math.sin(2*math.pi*i/55)
    rows.append([
        round(lam,       4),
        round(state[0],  6),
        round(state[1]+eps_r,  6),
        round(state[2]+eps_th, 6),
        round(state[3] % (2*math.pi), 6),
    ])
    state = rk4_step(kerr_rhs, lam, state, h, E, L, M, a)

# ── Write CSVs to current directory ──────────────────────────────────
out_dir = os.path.dirname(os.path.abspath(__file__))

coord_names = ['t_coord', 'r_coord', 'theta_coord', 'phi_coord']
for j, name in enumerate(coord_names):
    path = os.path.join(out_dir, f"isco_{name}.csv")
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([f"# Kerr ISCO geodesic a/M={a}"])
        w.writerow(['lambda', name])
        for row in rows:
            w.writerow([row[0], row[j+1]])
    print(f"Written {path}  ({N} samples)")

path = os.path.join(out_dir, "isco_all.csv")
with open(path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow([f"# Kerr ISCO geodesic a/M={a} — all four coordinates"])
    w.writerow(['lambda', 't_coord', 'r_coord', 'theta_coord', 'phi_coord'])
    for row in rows:
        w.writerow(row)
print(f"Written {path}  ({N} samples, all 4 coords)")

print(f"\nFirst 5 rows:")
print("lambda      t           r           theta       phi")
for row in rows[:5]:
    print("  ".join(f"{v:10.4f}" for v in row))
