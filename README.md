# kerr-black-hole

Analysis of Kerr black hole geodesic trajectories using signal geometry tools, with an X3D visualisation.

Each Boyer-Lindquist coordinate — t(λ), r(λ), θ(λ), φ(λ) — is treated as a scalar time series and fed through two companion tools: the [AVS Ground Tool](https://github.com/subunits/avsp) for augmented vector space analysis, and the [AVS Kähler Extension](https://github.com/subunits/kahler-ts) for Hermitian geometry and persistent homology.

---

## Black hole parameters

| Parameter | Value |
|---|---|
| Spin a/M | 0.7 |
| ISCO radius | 3.3931 M |
| Geodesic energy E | 0.8964 |
| Angular momentum L | 2.5865 |
| Samples | 300 (affine step h = 0.5) |
| Units | G = c = M = 1 |

---

## Files

    kerr_isco.py            RK4 Kerr geodesic integrator — generates trajectory CSVs
    isco_all.csv            All four coordinates combined (300 samples)
    isco_r_coord.csv        Radial coordinate r(λ) — input to Augmented Vector Space Ground Tool
    isco_t_coord.csv        Coordinate time t(λ) — input to Augmented Vector Space Ground Tool
    isco_theta_coord.csv    Polar angle θ(λ) — input to Augmented Vector Space Ground Tool
    isco_phi_coord.csv      Azimuthal angle φ(λ) — input to Augmented Vector Space Ground Tool
    avsp/Main.hs            Augmented Vector Space Ground Tool — kNN, OLS, anomaly detection
    kahler_isco.hs          AVS Kähler Extension — ISCO signal corpus
    black_hole_soft.x3d     Pure X3D 3.3 scene — softened Kerr black hole (a/M = 0.7)
    README.md               This file

---

## Requirements

- Python 3 — stdlib only (`math`, `csv`)
- GHC 9.x — `base` only, no Cabal or Stack required

---

## Step 1 — Generate trajectory data

    python3 kerr_isco.py

Writes to the current directory:

    isco_t_coord.csv        Coordinate time t(λ)
    isco_r_coord.csv        Radial coordinate r(λ)
    isco_theta_coord.csv    Polar angle θ(λ)
    isco_phi_coord.csv      Azimuthal angle φ(λ)
    isco_all.csv            All four coordinates combined

---

## Step 2 — AVS Ground Tool

Clone [avsp](https://github.com/subunits/avsp) and compile:

    ghc -O Main.hs -o avsp

Feed each coordinate CSV into the ground tool:

    LANG=C.UTF-8 ./avsp isco_r_coord.csv
    LANG=C.UTF-8 ./avsp isco_t_coord.csv
    LANG=C.UTF-8 ./avsp isco_theta_coord.csv
    LANG=C.UTF-8 ./avsp isco_phi_coord.csv

Each run produces kNN search, OLS next-step regression, Shannon entropy per augmented dimension, anomaly detection, and pairwise distance matrices (L2 and cosine) for that coordinate signal. Output also available on [play.haskell.org](https://play.haskell.org) — paste `Main.hs` and run.

---

## Step 3 — Kähler Extension

Compile and run `kahler_isco.hs` directly from this repo:

    ghc -O kahler_isco.hs -o kahler_isco
    LANG=C.UTF-8 ./kahler_isco

All four ISCO coordinates are embedded in `allSignals`. The primary signal `xs` is ISCO-r. Output covers Kähler condition verification, discrete curvature dω, Chern proxy comparison, symplectic form matrix, holomorphic kNN, and Vietoris-Rips persistent homology. Also runs on [play.haskell.org](https://play.haskell.org).

---

## Key results

### Chern proxy — symplectic integral Σ ω(pᵢ, pⱼ)

| Signal | Chern proxy | vs sine (1037.38) | Physical meaning |
|---|---|---|---|
| sine (ref) | 1037.38 | 100% | Pure periodic baseline |
| ISCO-θ | 987.57 | 95% | Near-sinusoidal polar oscillation, monotone phase |
| ISCO-r | 934.59 | 90% | Oscillatory radial perturbation, monotone accumulation |
| ISCO-φ | 512.45 | 49% | Azimuthal advance; 2π wrapping causes partial cancellation |
| lorenz | 436.26 | 42% | Chaotic Lorenz attractor — partial ω cancellation |
| ISCO-t | 23.56 | 2% | Linear coordinate time; near-collinear lag vectors |

### Curvature dω

All four ISCO coordinates show constant dω ≈ 0.0325 across all triples — uniform to three significant figures. The embedding is uniformly curved (not Lagrangian), consistent with a circular geodesic at fixed r_ISCO in a curved spacetime. This is the Kähler signature of orbital regularity.

### Symplectic matrix (ISCO-r, first 4 points)

    ω(u,v)    t=3      t=4      t=5      t=6
    t=3       0.000   +0.440   +0.870   +1.278
    t=4      -0.440    0.000   +0.440   +0.870
    t=5      -0.870   -0.440    0.000   +0.440
    t=6      -1.278   -0.870   -0.440    0.000

Uniform increment Δω = 0.440 per step — the orbital frequency of the ISCO encoded as symplectic phase advance in C^k.

### Persistent homology

H0 drops 15 → 1 as ε grows. H1 = 0 throughout. Entropy rises monotonically 0.0 → 3.25 bits.

---

## X3D scene

Open `black_hole_soft.x3d` in [X_ITE](https://create3000.github.io/x_ite/) (drag and drop) or Instant Reality. Cycle viewpoints with PageUp / PageDown — overview, equatorial, polar. Pure X3D 3.3, no external dependencies.

---

## Extending

Edit `kerr_isco.py` to analyse other orbital regimes:

- Plunging orbit: set `r < r_isco`, remove the `eps_r` stabilisation term
- Eccentric orbit: initialise with non-circular E and L values
- Different spin: change `a = 0.7` — r_ISCO recomputes automatically
- Photon sphere: set `r = 1.5 * M * (1 + sqrt(1 - a**2/M**2))` (approximate)

---

## Companion repositories

- [subunits/avsp](https://github.com/subunits/avsp) — AVS Ground Tool v3. CSV ingestion, kNN, OLS, anomaly detection. NASA ancillary software, EAR99.
- [subunits/kahler-ts](https://github.com/subunits/kahler-ts) — AVS Kähler Extension v2. Hermitian geometry, Chern proxy, persistent homology. NASA ancillary software, EAR99, TRL 2.

---

## Licence

Public domain — Unlicense. No warranty expressed or implied.
