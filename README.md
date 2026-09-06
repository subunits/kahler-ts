# AVS Kähler — Signal Manifold Analysis

## NASA Ancillary Software Compliance

**Software Type:** Ancillary — Research
**Classification:** Open Source
**Export Control:** EAR99 — No export restrictions apply
**Safety Critical:** No
**TRL:** 2 (Technology concept and/or application formulated)
**Language:** Haskell (GHC 9.x)
**Dependencies:** `base` only (no third-party packages)
**Platform:** Any platform with GHC; tested via [play.haskell.org](https://play.haskell.org)
**Companion project:** avs (Augmented Vector Space Ground Tool)

## Abstract

This software extends the Augmented Vector Space (AVS) framework into the domain of complex differential geometry. Each augmented time series vector in R^n is complexified by pairing consecutive lag dimensions into coordinates in C^(n/2). The resulting complex vector space admits a Hermitian inner product H = g + iω, where g is the Riemannian metric (recovering the original L2 distance) and ω is a symplectic 2-form encoding the phase relationships between lag dimensions. The software verifies the Kähler conditions on the complexified corpus, measures the curvature of the embedded signal manifold via discrete evaluation of dω, computes a Chern number proxy across all six built-in signal generators, and runs a Vietoris-Rips persistent homology filtration over the Hermitian distance to extract Betti numbers H0 and H1. The tool is intended for research into the topological and geometric structure of time series data and requires familiarity with Hermitian geometry and algebraic topology.

## Version History

| Version | Date | Author | Description |
|---------|------|--------|-------------|
| 1.0.0 | 2026-09-05 | (see Point of Contact) | Initial release: complexification, Kähler checks, curvature, Chern proxy, persistent homology |

## 1. Mathematical Background

### 1.1 Prerequisites

The reader is assumed to be familiar with the following:

- Real and complex inner product spaces
- Riemannian geometry (metric tensor, geodesics)
- Symplectic geometry (symplectic form, Lagrangian submanifolds)
- Kähler manifolds (Hermitian manifolds with closed symplectic form)
- Basic algebraic topology (homology groups, Betti numbers)
- Persistent homology and the Vietoris-Rips filtration

### 1.2 From AVS to Kähler

The AVS ground tool constructs an augmented vector v(t) in R^n for each time step t by concatenating lag features, finite differences, and rolling statistics. This is a flat real vector space with the standard L2 metric.

The Kähler extension proceeds in three steps.

First, complexification. Consecutive lag pairs (lag_{2k}, lag_{2k+1}) are identified with the real and imaginary parts of a complex coordinate z_k = lag_{2k} + i·lag_{2k+1}. Derivative and rolling statistic dimensions are embedded as purely real complex numbers. The result is a point in C^m where m = ceil(p/2) + d + 2, and p, d are the lag window and derivative order respectively.

Second, the Hermitian inner product. The standard Hermitian form H(u,v) = Σ conj(u_i)·v_i decomposes as H = g + iω where g(u,v) = Re H(u,v) is the Riemannian metric and ω(u,v) = Im H(u,v) is the symplectic 2-form. Since Re H(v,v) = Σ |v_i|² = ‖v‖², the Riemannian metric exactly recovers the original L2 norm. The complexification preserves the metric.

Third, the complex structure. The map J: C^m → C^m defined by J(z) = iz acts as a 90-degree rotation in each complex plane. It satisfies J² = -I, g(Ju,Jv) = g(u,v), and the compatibility condition ω(u,v) = g(Ju,v). Together with the closedness of ω (dω = 0 on the flat C^m), these conditions define a Kähler manifold. The flat space C^m with the standard Hermitian metric is the simplest example of a Kähler manifold.

### 1.3 Embedded Signal Manifold

The augmented points {v(t)} trace out a curve (or low-dimensional manifold) embedded in C^m. This embedded manifold is generally not flat even though the ambient C^m is. The discrete exterior derivative dω evaluated on triples of consecutive signal points measures the local curvature of the embedding. Specifically, dω(p,q,r) = ω(q-p,r-p) + ω(r-q,p-q) + ω(p-r,q-r). For a flat embedding this vanishes; non-zero values quantify how far the signal manifold deviates from a Lagrangian submanifold of C^m.

### 1.4 Chern Number Proxy

The Chern number of a complex vector bundle is a topological invariant computed as the integral of the curvature 2-form over the base manifold. The discrete proxy computed here, Σ_{i<j} ω(p_i,p_j), approximates this integral over the corpus of augmented points. For signals with monotone phase progression (sine wave) this sum is large because ω accumulates without cancellation. For chaotic signals whose trajectories cross themselves in lag space (logistic map, random walk) partial cancellation reduces the net integral.

### 1.5 Persistent Homology

The Vietoris-Rips filtration builds a simplicial complex over the corpus at each distance threshold ε. Vertices are augmented points; edges connect pairs within distance ε; triangles fill triples within mutual distance ε. The Betti numbers β0 (connected components) and β1 (independent 1-cycles) track the topology of this complex as ε grows.

The homology corpus consists of two guaranteed-separated clusters: 8 points from a sine wave centred at 0 and 8 points from the same sine wave DC-shifted by +10. The DC shift moves the lag vectors to a completely separate region of C^k — the inter-cluster gap is approximately 10√dim, which exceeds the maximum intra-cluster distance by construction. This guarantees β0 starts at 2 and drops to 1 cleanly as ε crosses the inter-cluster gap. Filtration steps follow sorted pairwise distances rather than uniform spacing so every topology-changing event is captured. Each step also reports the Shannon entropy of the active edge distance distribution, which rises as the complex gains structurally diverse edges. Connected components are computed via a stateless union-find over plain lists, avoiding state-threading issues while remaining correct for the corpus sizes used.

## 2. Software Description

### 2.1 Architecture

The program is a single Haskell module (Main) with the following layers. Signal generators produce raw sample sequences. The AVS augmentation layer (inherited from the companion project) constructs real augmented vectors. The complexification layer lifts these into C^m. The geometry layer computes the Hermitian metric, symplectic form, complex structure J, and curvature. The topology layer runs the Vietoris-Rips filtration and extracts Betti numbers. All results are reported to stdout.

### 2.2 Key Data Types

| Type | Description |
|------|-------------|
| `C` | Complex number: `C { re :: Double, im :: Double }` |
| `CVec` | `[C]` — a complex feature vector in C^m |
| `AugConfig` | Augmentation hyperparameters (inherited from AVS) |
| `AugPoint` | Real augmented point: time index + real vector |
| `KahlerPoint` | Kähler-lifted point: time index + real vector + complex vector |
| `Parents` | `[Int]` — stateless union-find parent array for VR connected components |

### 2.3 Algorithms

#### Complexification

Consecutive lag pairs are identified as complex coordinates. Purely real dimensions (diffs, rolling stats) are embedded as C x 0. Time complexity O(n·m) where n is the number of augmented points and m is the complex dimension.

#### Hermitian inner product

H(u,v) = Σ conj(u_i)·v_i computed in O(m) per pair.

#### Kähler condition verification

Four conditions checked per consecutive pair: J²=-I, g(Ju,Jv)=g(u,v), ω skew-symmetric, ω(u,v)=g(Ju,v). All are O(m) per pair.

#### Curvature (discrete dω)

Evaluated on consecutive triples. O(m) per triple.

#### Chern proxy

O(n²·m) over the full corpus.

#### Vietoris-Rips persistent homology

Edges O(n²), triangles O(n³). Restricted to the first 12 points by default to keep runtime tractable on play.haskell.org. For larger corpora, compile locally with ghc -O.

## 3. Inputs and Outputs

### 3.1 Inputs

All parameters are set by editing constants in main.

| Parameter | Variable | Default | Description |
|-----------|----------|---------|-------------|
| Lag window | `lagWindow` | 4 | Number of lag dimensions (even recommended for clean pairing) |
| Derivative order | `derivOrder` | 1 | 0, 1, or 2 |
| Rolling window | `rollingWin` | 4 | Window for μ and σ features |
| Signal | `xs` | `lorenzWave 60` | Input time series |
| Homology points | `kPts12` | first 12 | Points used in VR filtration (increase carefully) |

Built-in signals:

| Function | Type | Notes |
|----------|------|-------|
| `sineWave n` | Periodic | Near-zero Chern proxy |
| `logisticMap n` | Chaotic | Large Chern proxy, high entropy |
| `lorenzWave n` | Chaotic | Large curvature flux; intermediate Chern proxy |
| `sawtoothWave n` | Quasi-periodic | Intermediate curvature |
| `stepWave n` | Piecewise+noise | Step discontinuities visible in ω |
| `randomWalk n` | Non-stationary | Growing norm, diffuse ω |

### 3.2 Outputs

All output to stdout. No files are written (unlike the AVS ground tool). Sections produced:

| Section | Description |
|---------|-------------|
| Kähler condition verification | J²=-I, metric preservation, skew-symmetry, compatibility |
| Embedded manifold curvature | dω per triple; magnitude = local curvature flux |
| Chern proxy comparison | All six signals ranked by discrete ω integral |
| Symplectic form matrix | ω(p_i, p_j) for first 6 points |
| Holomorphic 5-NN | Nearest neighbours in Hermitian distance |
| Norm comparison | Real L2 vs Hermitian norm (ratio = 1.0 confirms correctness) |
| Persistent homology | H0 and H1 Betti numbers across VR filtration |

## 4. Usage

### 4.1 Online

Open https://play.haskell.org, paste Main.hs, press Run. The homology corpus is 16 points (8 per cluster) and the filtration runs in O(n³) triangle time — fast enough to complete within the playground timeout.

### 4.2 Local

```bash
ghc -O Main.hs -o avs-kahler
./avs-kahler
```

The homology corpus size is controlled by the `stride` function in main. Increasing the number of points raises runtime as O(n³).

### 4.3 Interpreting the output

The Kähler checks will always show all ✓ — the flat C^m space always satisfies these conditions. The interesting outputs are the curvature (dω on triples, non-zero for chaotic signals), the Chern proxy comparison across all six signals (sine largest, random walk smallest), and the persistent homology filtration showing H0 dropping from 15 to 1 as the two clusters merge, with entropy rising monotonically from 0.0 to ~3.25 bits.

## 5. Relationship to AVS Ground Tool

This project is the research extension of the AVS Ground Tool (avs). The ground tool is suitable for operational use against telemetry data: it ingests CSV files, detects anomalies, and exports results. This project is not intended for operational use. It is a mathematical workbench for studying the geometry and topology of the signal manifold.

The two projects share the augmentation layer (AugConfig, augment, signal generators) but diverge at the point where the AVS ground tool applies kNN and OLS and this project applies Hermitian geometry and persistent homology.

## 6. Limitations

The Vietoris-Rips filtration is O(n³) in the number of triangles and is restricted to small corpora. The Chern proxy is a discrete approximation to a continuous integral and is sensitive to corpus size and signal length. The persistent homology implementation uses a simplified Euler-characteristic H1 estimator rather than a full boundary matrix reduction; it is suitable for exploratory analysis but not for publication-quality topological data analysis. The complexification pairs consecutive lag dimensions, which is natural but not unique; other pairings (e.g. by frequency content via DFT) would produce different complex structures. All checks are performed on the flat ambient C^m; the intrinsic geometry of the embedded signal manifold is only probed indirectly through dω and the VR filtration.

## 7. Point of Contact

| Field | Value |
|-------|-------|
| Author | — |
| Organisation | — |
| Email | — |
| Distribution | Unlimited (EAR99) |

This README was prepared in accordance with NASA NPR 2210.1C requirements for ancillary software. The software has not been subjected to NASA IV&V and is not intended for flight or safety-critical use. TRL 2 designation reflects that the mathematical concept has been formulated and demonstrated in software but has not been validated against mission data.

## 8. Licence

This software is released into the public domain under the Unlicense (https://unlicense.org). No warranty is expressed or implied.
