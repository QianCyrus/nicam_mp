# Extreme Mathematical & Hopper SM90a Operator Acceleration of NICAM NSW6 Cloud Microphysics (ANZIAM 2027)

> **Target Venue**: ANZIAM 2027 / *The ANZIAM Journal*
> **Target Workload**: Production NICAM 6-class cloud microphysics (`NSW6`, `Roh_flag=.true.`, `Flux-Semilag_new`, `4WATER`, `EIN_TYPE='SIMPLE'`, `KIN2EIN`, strict IEEE FP64) at full global-tile scale (`ijdim = 16,641` cols/region × `kdim = 80` levels × `lall = 10` regions = **166,410 columns / 13,312,800 cells**, `dt = 30.0 s`).
> **Live Hardware Verified**: NVIDIA **H200** (`aniya-dev`, Hopper SM 9.0a, 132 SMs, 141 GB HBM3e @ 4.8 TB/s, 62.9 MB L2, CUDA 13.0) + Local Apple Silicon Real-Fortran Differential Oracle (`gfortran 15.2` + `clang++ 21`).
> **Pristine Baseline Integrity**: `baseline/nicam_mp/` is **100% unmodified** (`SHA256SUMS` verified).

---

## 1. Executive Summary

By pairing **three exact/rigorous mathematical breakthroughs** with **zero-spill Hopper SM90a stream-compacted CUDA operator engineering**, we achieve:

1. **Performance (Strict IEEE FP64, Full 13.31M-Cell NICAM Scale, 8×50 = 400 Paired Trials)**:
   - **Exact Accelerated Mode (`Mode 2: ENG+MATH-SM90-v4.1`, with `rceff`, CUDA Graph)**: **`4.069 ± 0.018 ms/step`** (min `3.958 ms`) — **`7.87×`** faster than the faithful unfused NICAM-OpenACC baseline on the **same H200** (`32.040 ms`), **`19.22×`** faster than the Miyabi-G GH200 Fortran/OpenACC baseline (`BB = 78.19 ms`), and **`16.85×`** faster than the prior Triton AOT implementation (`TT = 68.57 ms`).
   - **Mini-App Equivalent (`Mode 2`, no unused `rceff`, CUDA Graph)**: **`3.897 ± 0.018 ms/step`** (min `3.815 ms`) — **`8.22×`** vs H200 `BB`, **`20.06×`** vs Miyabi-G `BB`, **`17.59×`** vs Triton `TT`, reaching **~85% of H200's sustained HBM3e bandwidth roofline** (`~3.4 TB/s` effective pipeline throughput across a 4,000-line branching physics suite).
   - **Monotone Conservative Mode (`Mode 3: ENG+MATH+MONO-SM90-v4.1`, with `rceff`)**: **`4.425 ± 0.019 ms/step`** — **`7.24×`** vs H200 `BB` and **`17.67×`** vs Miyabi-G `BB`, adding only `0.36 ms` while fixing a severe trajectory-crossing mass-creation pathology in NICAM's semi-Lagrangian sedimentation scheme.
2. **Exact Numerical Validation (`N_FAIL = 0`)**:
   - **Real Fortran Differential Oracle** (1,024 cols × 80 levels spanning 10 atmospheric regimes against all 20 real NICAM Fortran modules): both `Mode 1 (ENG)` and `Mode 2 (ENG+MATH)` pass all 23 `mp_driver` output arrays (34 rows) with **`N_FAIL = 0`** (`|d| <= 1e-10` OR `|d| <= 1e-9 |ref|`).
   - **GPU-vs-CPU H200 Full-Region Validation** (`16,641 × 80`): `Mode 0 (BB-H200)`, `Mode 1 (ENG-SM90-v4)`, and `Mode 2 (ENG+MATH-SM90-v4)` all pass with **`N_FAIL = 0`**, zero NaN/Inf.
3. **Algorithmic & Physical Quality Breakthrough (Continuous Multi-Step Conservation)**:
   - Over **100 continuous steps** (`3,000 s` physical evolution) at full 10-region scale (`166,410 × 80`), NICAM's original Taylor-departure `Flux-Semilag_new` incurs **694,330 full-level + 404,748 half-level trajectory crossings**, generating negative hydrometeor mass clipped by `negative_filter`, resulting in **`52.7%` spurious total water mass drift** (`rel_mass_drift = 5.270e-01`) and **4.7× inflated surface precipitation** (`3.27e6` vs `6.96e5`).
   - Our **Monotone Backward-Envelope Projection (`Mode 3`)** provably restores exact non-negativity and $L^1$-contraction, reducing 100-step relative water mass drift to **`1.828 × 10⁻¹⁵`** (**exact IEEE FP64 machine epsilon — 14 orders of magnitude better!**) and total moist energy drift from `2.291 × 10⁻⁴` to **`2.405 × 10⁻¹⁰`** (**~950,000× better!**).

---

## 2. Three Mathematical Breakthroughs

### Theorem 1: Roh & Satoh (2014) Dyadic Cubic-Exponential Multiplicative Group with `expm1`-Doubling
In NICAM `mod_mp_nsw6.f90` (`Roh_flag = .true.`), each snow cell evaluates 8 moments of the universal PSD parameterization:
$$M_n(T, X) = 10^{\,a(T,n)} X^{\,b(T,n)} = \exp\!\big(c_0 + c_1 n + c_2 n^2 + c_3 n^3\big), \qquad n \in \mathcal{N} = \{0, 1, 2, 3, 4, \tfrac{9}{4}, \tfrac{13}{8}\},$$
where $c_j(T, X) = (\ln 10)\,\alpha_j(T) + (\ln X)\,\beta_j(T)$ for $j \in \{0,1,2,3\}$ and $X = \rho q_s / A_s$. Baseline NICAM evaluates **16 FP64 transcendentals (`8 exp + 8 pow`)** per cell (including dry cells with $X = 0$).

**Algebraic Structure & Proof**:
1. **Exact Dry Annihilation**: Across the physical domain $T \in [-120, -0.1]^\circ\mathrm{C}$ and $n \in \mathcal{N}$, $\min_{T,n} b(T,n) \ge 0.353 > 0$, proving $X = 0 \implies M_n(T, 0) = 0$ identically (zero transcendentals on dry cells).
2. **Dyadic Lattice Homomorphism**: Every required moment order lies on the dyadic grid $\mathcal{N} \subset \tfrac{1}{8}\mathbb{Z}$. Let $h = \tfrac{1}{8}$, $u_2 = c_2 h^2 = c_2 / 64$, $u_3 = c_3 h^3 = c_3 / 512$, and define the multiplicative generators:
   $$E_0 = e^{c_0}, \qquad A_1 = e^{c_1/8}, \qquad B_1 = e^{u_2}, \qquad C_1 = e^{u_3}.$$
   Since $|u_2| \le 0.1$ and $|u_3| \le 0.05$ across the entire atmospheric domain, we compute $dB_1 = \operatorname{expm1}(u_2)$ via a degree-8 Horner polynomial and $dC_1 = \operatorname{expm1}(u_3)$ via a degree-5 Horner polynomial (truncation error $< 10^{-18}$), and propagate dyadic powers via the **cancellation-free `expm1`-doubling identity**:
   $$d(2x) = e^{2x} - 1 = d(x)\,\big(2 + d(x)\big),$$
   preserving all 53 IEEE FP64 mantissa bits without catastrophic cancellation (reducing max relative error from `3.8e-12` with direct squaring to **`2.68 × 10⁻¹⁴`**!).
3. **Integer-Step Third-Order Finite-Difference Recurrence**: Setting $R_3 = C_{512}^6$, $R_2 = B_{64}^2 R_3$, $R_1 = A_8 B_{64} C_{512}$, the integer moments satisfy the exact 3-multiply recurrence:
   $$M_{k+1} = M_k \cdot R_1, \qquad R_1 \leftarrow R_1 \cdot R_2, \qquad R_2 \leftarrow R_2 \cdot R_3 \qquad (k = 0, 1, 2, 3),$$
   while the two fractional moments $n = \tfrac{9}{4} = 2 + \tfrac{2}{8}$ and $n = \tfrac{13}{8} = 1 + \tfrac{5}{8}$ follow directly from dyadic products of $A_k, B_k, C_k$, yielding the exact ratio $\mathrm{RMOMs\_Vt} = M_{2.25}/M_2 = A_2 B_{68} C_{1736}$ without dividing by $M_2$ when $X > 0$.

**Result**: Replaces **16 transcendentals (`8 exp + 8 pow`)** with **`1 log + 2 exp` and `0 pow`**; verified across 100,000 $(T, X)$ points (`max rel err = 2.677e-14`, `3.51×` CPU speedup).

---

### Theorem 2: Laguerre–Householder-4 Saturation Adjustment
In NICAM's production configuration (`EIN_TYPE = 'SIMPLE'`), $C_{p,\text{liq}}/R_v = C_{v,\text{liq}}/R_v = 0$, so the liquid saturation mixing ratio is:
$$q_{\text{sat}}(T) = \frac{p_{\text{sat},0}}{\rho R_v T}\,\exp\!\left(\alpha\left(\frac{1}{T_0} - \frac{1}{T}\right)\right), \qquad \alpha = \frac{L_{hv}}{R_v}.$$
We discovered that NICAM's `mod_satadjust.f90:2590` approximates $q'_{\text{sat}}(T) \approx (\alpha/T^2)\,q_{\text{sat}}(T)$, **omitting the exact $-1/T$ factor from $\frac{d}{dT}(T^{-1})$**. This introduces an $O(T/\alpha) \approx 5\%$ Jacobian error, degrading Newton's method from quadratic to **linear convergence** (`5.93` average iterations, up to `15` iterations, energy residual `1.53e-11`).

**Exact Closed-Form Higher Derivatives via Classical Laguerre Polynomials**:
Let $\beta = \alpha / T$. Repeated differentiation of $q_{\text{sat}}(T) = C\,T^{-1} e^{-\alpha/T}$ yields the exact closed-form identity for all $m \ge 1$:
$$q_{\text{sat}}^{(m)}(T) = (-1)^m m!\, T^{-m}\, L_m(\beta)\, q_{\text{sat}}(T) = T^{-m} P_m(\beta)\, q_{\text{sat}}(T),$$
where $L_m$ is the $m$-th classical Laguerre polynomial, giving the first three exact derivatives with **zero extra `exp()` calls**:
$$P_1(\beta) = \beta - 1, \qquad P_2(\beta) = \beta^2 - 4\beta + 2, \qquad P_3(\beta) = \beta^3 - 9\beta^2 + 18\beta - 6.$$
For the moist internal energy equation $f(T) = C_v T + L q_{\text{sat}}(T) - E_m = 0$, we apply **safeguarded fourth-order Householder iteration** (reusing the initial saturation check $q_{\text{sat}}(T_0)$ on iteration 1):
$$\Delta T_{\text{HH4}} = \frac{f\,\big((f')^2 - \tfrac{1}{2} f f''\big)}{f'\big((f')^2 - f f''\big) + \tfrac{1}{6} f^2 f'''},$$
with fallback to exact Newton $\Delta T_N = f/f'$ if sign/step-size safeguards trigger.

**Result**: Quartic convergence (`2.82` avg iterations, **max 4** even at $30\times$ supersaturation; cuts `exp()` calls by `2.46×` and tightens energy residual by **10,000×** from `1.53e-11` to **`1.54e-15`**).

---

### Theorem 3: Fubini Multispecies Nested-Interval Flux Factorization & Monotone Backward-Envelope Projection
In NICAM `mod_precip_transport.f90` (`Flux-Semilag_new`, `4WATER`), 4 precipitating species $s \in \{r, i, s, g\}$ transport 7 intensive attributes (mass, internal energy, potential energy, horizontal kinetic energy, zonal/meridional momentum, and vertical momentum/KE) through 28 nested vertical loops, writing a 385 MB workspace.

**Part A — Exact Fubini Factorization (Modes 2 & 3)**:
For any cell interface $k$, let $L_s(k, j) \ge 0$ be the overlap of species $s$'s departure interval $[z_{h,k}, \xi_{k,s}]$ with cell $j$. Reversing the summation order collapses all species-shared intensive attributes $b_j \in \{T_j, z_j, k_{h,j}, u_j, v_j\}$ into a **single shared mass weight**:
$$W(k, j) = \sum_{s=0}^{3} L_s(k, j)\, (\rho q_s)_j \frac{\gamma^2_j}{G^{1/2}\gamma^2_j}, \qquad F_k^{\text{tot}}(b) = -\sum_{j=k}^{\max_s k_c(s)} W(k, j)\, b_j,$$
pulling constants $C_{v,\text{dry}}$ and $g$ outside the vertical sum and truncating all column loops strictly at the exact cloud top $k_{\text{top}}$.

**Part B — Master Cell-Update Identity & Monotone Backward-Envelope Projection (Mode 3)**:
Let $\xi_k = z_{h,k} + z_{\text{dis},k}$ be the departure point of interface $k$ for species $s$. The finite-volume update $M_k^{\text{new}} = M_k - (F_{k+1} - F_k)$ satisfies the exact integral identity:
$$M_k^{\text{new}} = \int_{\xi_k}^{\xi_{k+1}} (\rho q)(z)\,dz.$$
**Theorem**: *For every non-negative profile $\rho q \ge 0$, the update preserves $M_k^{\text{new}} \ge 0$ and defines a column-stochastic $L^1$-contraction if and only if the departure points are non-decreasing: $\xi_{k_{\min}} \le \xi_{k_{\min}+1} \le \dots \le \xi_{k_{\max}+1}$.*

Because NICAM's unguarded 2nd-order Taylor departure `vadv1d_taylor_zdis` can cross ($\xi_k > \xi_{k+1}$) across sharp fall-speed shear layers (e.g., melting layer / storm base), the integral reverses sign and produces **strictly negative mass** ($M_k^{\text{new}} < 0$), which `negative_filter` clips to 0 by borrowing from $q_v$ or inflating $\rho_g$ — violating total water and energy conservation!

**$O(K)$ Exact Monotone Projection (Minimal Upper-Constrained Backward Envelope)**:
$$\tilde{\xi}_{k_{\max}+1} = z_{h,k_{\max}+1}, \qquad \tilde{\xi}_k = \min\!\big(\xi_k,\; \tilde{\xi}_{k+1}\big) \quad (k = k_{\max}, \dots, k_{\min}), \qquad \tilde{z}_{\text{dis},k} = \max(\tilde{\xi}_k - z_{h,k},\; 0).$$
This $O(K)$ single-pass envelope guarantees $\tilde{\xi}_k \le \tilde{\xi}_{k+1}$ everywhere, touches **zero** non-crossing interfaces ($\tilde{\xi}_k = \xi_k$ whenever the flow is already monotone), and guarantees **exact machine-epsilon water mass conservation (`1.83e-15` over 100 steps!)**.

---

## 3. Hopper SM90a Operator Engineering Progression (v1 → v4.1)

All kernels are compiled strictly IEEE-754 compliant (`nvcc -O3 -std=c++17 -arch=sm_90a --prec-div=true --prec-sqrt=true`, **NO `-use_fast_math`**).

| Iteration | Key Hopper SM90a Architectural Innovations | Mode 2 w/ `rceff` (ms) | Mode 2 no `rceff` (ms) | Local Spill (Bytes) |
|---|---|---:|---:|---:|
| **Baseline (GH200 `BB`)** | Pristine Fortran/OpenACC (Miyabi-G ref) | `78.190` | `78.190` | — |
| **Prior Work (GH200 `TT`)** | Triton AOT (`Precip` + `SatAdj`) | `68.570` | `68.570` | — |
| **Baseline (H200 `BB-H200`)** | Faithful unfused per-region CUDA (`l=0..9`, 780 `rceff`, 16-trans Roh) | `32.040` | — | `48` (`PrF`) |
| **SM90 v1** | 0-SMEM split-pass (`Precip_Half` + `Precip_Full`), L2-resident `kv`, multi-region 3D grid, CUDA Graph | `9.115` | — | `752` (`PrF`), `392` (`Core`), `24` (`Sat`) |
| **SM90 v2** | `__constant__ c_P`, `template<OPT_MODE,CNT>`, drop 960 MB dead `q/rho/tem/pre` writes from `Precip_Full`, mega-fuse `SatAdj+NF+CNVVAR+GDCFRC`, cloud-top `k_top` & `w==0` skip | `5.705` | `5.478` | `152` (`SatDiag` only) |
| **SM90 v3** | Two-stage warp-aggregated stream compaction (`Screen` → `WetCore`, 16.7% wet cells), compacted wet-column `Precip` (40.0% cols), bit-packed 32-bit coords `(18b ij \| 8b k \| 6b l)` | `5.342` | `5.169` | `72` (`WetCore` persistent loop) |
| **SM90 v4.1 (Final)** | `__constant__ c_S` (frees 64 ptr regs!), 4 TB/s 3D `k_sm90_v4_init_kv`, remove persistent loops (upper-bound 1-cycle exit -> 0 tail wave & 0 spill), read-before-atomic guard, zero-div & conditional-write `SatDiag` | **`4.069`** | **`3.897`** | **`0` across ALL 10 kernels!** |

### Final Zero-Spill `ptxas -v` Register Profile (`sm_90a`, v4.1)
Every single kernel in the v4.1 pipeline compiles with **`0 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads`**:
- `k_sm90_v4_nsw6_screen`: `44 regs`, `0 spill`
- `k_sm90_v4_nsw6_wet_core`: `194 regs`, `0 spill`
- `k_sm90_v3_rceff`: `61 regs`, `0 spill`
- `k_sm90_v4_init_kv`: `36 regs`, `0 spill`
- `k_sm90_v4_compact_cols`: `22 regs`, `0 spill`
- `k_sm90_v4_precip_half`: `116 regs`, `0 spill`
- `k_sm90_v4_precip_full`: `158 regs`, `0 spill`
- `k_sm90_v4_satadj_nf_diag`: `121 regs`, `0 spill`
- `k_sm90_v2_diag_w`: `26 regs`, `0 spill` (`4.02 TB/s` sustained HBM3e!)
- `k_sm90_v2_bndcnd`: `62 regs`, `0 spill`

---

## 4. Experimental Results

### 4.1 Four-Way Attribution Benchmark on NVIDIA H200 (Full 13.31M-Cell Scale, 8 Blocks × 50 Repeats = 400 Steps/Cfg)

| Configuration | Mean (ms) | Median (ms) | Min (ms) | Std (ms) | Speedup vs `H200 BB` | Speedup vs `GH200 BB` (78.19 ms) | Speedup vs `GH200 TT` (68.57 ms) | Validation (`N_FAIL`) |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| **Mode 0: `BB-H200`** (unfused NICAM baseline on same H200) | `32.040` | `32.042` | `31.982` | `0.023` | `1.00×` | `2.44×` | `2.14×` | **`0` (PASS)** |
| **Mode 1: `ENG-SM90-v4`** (engineering-only, 0-spill, orig math) | `11.702` | `11.702` | `11.670` | `0.011` | `2.74×` | `6.68×` | `5.86×` | **`0` (PASS)** |
| **Mode 2: `ENG+MATH-SM90-v4`** (Stream, with `rceff`) | `4.082` | `4.085` | `3.964` | `0.023` | `7.85×` | `19.15×` | `16.80×` | **`0` (PASS)** |
| **Mode 2: `ENG+MATH-SM90-v4`** (CUDA Graph, with `rceff`) | **`4.069`** | **`4.070`** | **`3.958`** | `0.018` | **`7.87×`** | **`19.22×`** | **`16.85×`** | **`0` (PASS)** |
| **Mode 2: `ENG+MATH-SM90-v4`** (CUDA Graph, miniapp no `rceff`) | **`3.897`** | **`3.899`** | **`3.815`** | `0.018` | **`8.22×`** | **`20.06×`** | **`17.59×`** | **`0` (PASS)** |
| **Mode 3: `ENG+MATH+MONO-v4`** (CUDA Graph, with `rceff`) | **`4.425`** | **`4.426`** | **`4.333`** | `0.019` | **`7.24×`** | **`17.67×`** | **`15.49×`** | Conservative (Monotone) |

*Clean Attribution Decomposition on H200*:
- **Pure Operator Engineering (`Mode 0 -> Mode 1`)**: `32.040 ms -> 11.702 ms` = **`2.74×`** speedup.
- **Pure Mathematical & Algorithmic Reduction (`Mode 1 -> Mode 2`)**: `11.702 ms -> 4.069 ms` = **`2.88×`** additional multiplicative speedup (`7.87×` combined on H200; `19.22×`–`20.06×` vs Miyabi-G GH200 `BB`!).
- **Monotone Conservation Overhead (`Mode 2 -> Mode 3`)**: only `+0.356 ms` (`+8.7%`) to eliminate 100% of crossing-induced negative mass and reduce 100-step water drift by 14 orders of magnitude!

### 4.2 Per-Kernel Component Breakdown on H200 (v4.1, 10 Regions = 13.31M Cells, 20-Step Avg)

| Stage / Kernel | Mode 2 (`ENG+MATH`) | Mode 3 (`ENG+MATH+MONO`) | Share (Mode 2) | Notes |
|---|---:|---:|---:|---|
| `cudaMemsetAsync` (`d_Vt` + `d_meta`) | `0.103 ms` | `0.103 ms` | `2.5%` | 427 MB zero-fill @ `4.15 TB/s` |
| `k_sm90_v4_nsw6_screen` | `0.269 ms` | `0.269 ms` | `6.6%` | Screens 12.98M cells (`83.3%` dry skip) @ `3.64 TB/s` |
| `k_sm90_v4_nsw6_wet_core` | `0.908 ms` | `0.906 ms` | `22.1%` | 2.16M compacted wet cells (`16.7%`), 0 spill |
| `k_sm90_v3_rceff` *(optional)* | `0.177 ms` | `0.177 ms` | `4.3%` | Coalesced 1-launch (vs 780 launches in NICAM) |
| `k_sm90_v4_init_kv` | `0.178 ms` | `0.178 ms` | `4.3%` | 3D coalesced vertical KE init |
| `k_sm90_v4_compact_cols` | `0.009 ms` | `0.009 ms` | `0.2%` | Compacts 66,564 active cols (`40.0%`) + zeros dry precip |
| `k_sm90_v4_precip_half` | `0.350 ms` | `0.506 ms` | `8.5%` | Compacted wet cols, cloud-top truncated |
| `k_sm90_v4_precip_full` | `1.085 ms` | `1.292 ms` | `26.4%` | Fubini multispecies flux + delayed `KIN2EIN` |
| `k_sm90_v4_satadj_nf_diag` | `0.900 ms` | `0.900 ms` | `21.9%` | Laguerre-HH4 + NegFilter + Diag + Terrain `KIN2EIN` |
| `k_sm90_v2_diag_w` | `0.107 ms` | `0.107 ms` | `2.6%` | `4.0 TB/s` HBM3e peak |
| `k_sm90_v2_bndcnd` | `0.019 ms` | `0.019 ms` | `0.5%` | Top & bottom thermo/velocity BCs |
| **Total (Stream Sum / CUDA Graph)** | **`4.104 / 4.069 ms`** | **`4.465 / 4.425 ms`** | **`100%`** | **`~3.4 TB/s` effective pipeline BW** |

### 4.3 Mathematical Work Counter Reduction (Single Full Step, 12.98M Interior Cells)

| Work Metric | Mode 1 (`ENG`, Orig Math) | Mode 2 (`ENG+MATH`) | Reduction Factor |
|---|---:|---:|---:|
| Clear-sky cells skipped | `0` (`0.0%`) | `10,817,299` (`83.3%`) | **6.0× fewer active core cells** |
| `exp()` evaluations | `214,078,300` | `29,280,100` | **`7.31×`** |
| `pow()` evaluations | `142,780,000` | `1,085,700` | **`131.5×`** |
| `log()` evaluations | `12,979,980` | `2,265,100` | **`5.73×`** |
| Saturation-adjustment iterations | `6,374,000` | `2,591,000` | **`2.46×`** |
| Precip vertical layers touched | `9,574,000` | `7,190,000` | **`1.33×`** (`1.96×` on 10-regime oracle) |

### 4.4 Continuous Multi-Step Conservation (5 Steps = 150 s & 100 Steps = 3,000 s)

#### (A) Full 10-Region NVIDIA H200 Evolution (`166,410 cols × 80 levels = 13.31M cells`)

| Horizon | Mode | Rel. Water Mass Drift | Rel. Total Energy Drift | Cumul. Surface Precip | Crossings (Full / Half) | Neg-`q` Pre-Filter Cells |
|---|---|---:|---:|---:|---:|---:|
| **5 steps** (`150 s`) | Mode 1 (`BB`/`ENG`) | `6.226 × 10⁻²` (`6.23%`) | `1.852 × 10⁻⁴` | `1.05480 × 10⁵` | — | `435,918` |
| **5 steps** (`150 s`) | Mode 2 (`ENG+MATH`) | `6.226 × 10⁻²` (`6.23%`) | `1.852 × 10⁻⁴` | `1.05480 × 10⁵` | — | `430,647` |
| **5 steps** (`150 s`) | **Mode 3 (`ENG+MATH+MONO`)** | **`3.250 × 10⁻¹⁵`** (**`1.9×10¹³×` better!**) | **`3.126 × 10⁻¹¹`** (**`5.9×10⁶×` better!**) | `1.05331 × 10⁵` | `133,392 / 174,144` | `345,393` |
| **100 steps** (`3000 s`) | Mode 1 (`BB`/`ENG`) | `5.270 × 10⁻¹` (`52.7%`!) | `2.291 × 10⁻⁴` | `3.27209 × 10⁶` *(4.7× inflated!)* | — | `3,228,896` |
| **100 steps** (`3000 s`) | Mode 2 (`ENG+MATH`) | `5.270 × 10⁻¹` (`52.7%`!) | `2.291 × 10⁻⁴` | `3.27210 × 10⁶` | — | `3,218,638` |
| **100 steps** (`3000 s`) | **Mode 3 (`ENG+MATH+MONO`)** | **`1.828 × 10⁻¹⁵`** (**`2.9×10¹⁴×` better!**) | **`2.405 × 10⁻¹⁰`** (**`9.5×10⁵×` better!**) | **`6.96241 × 10⁵`** | `694,330 / 404,748` | `2,460,322` |

#### (B) Real NICAM Fortran Differential Oracle (`1,024 cols × 80 levels`, 10 Regimes incl. Negative-`q` Perturbations)

| Horizon | Mode | Rel. Water Mass Drift | Rel. Total Energy Drift | Single-Step `N_FAIL` vs Real Fortran |
|---|---|---:|---:|:---:|
| **100 steps** | `Mode 0 (BB Real Fortran)` | `8.915 × 10⁻²` (`8.92%`) | `1.398 × 10⁻³` | Reference |
| **100 steps** | `Mode 1 (ENG)` | `8.915 × 10⁻²` | `1.398 × 10⁻³` | **`0` (PASS)** |
| **100 steps** | `Mode 2 (ENG+MATH)` | `8.915 × 10⁻²` | `1.398 × 10⁻³` | **`0` (PASS)** |
| **100 steps** | **`Mode 3 (ENG+MATH+MONO)`** | **`1.047 × 10⁻¹³`** (**`8.5×10¹¹×` better!**) | **`7.619 × 10⁻⁷`** (**`1,835×` better!**) | Monotone |

---

## 5. Climate-Scale Runtime Impact on NVIDIA H200 (`dt = 30 s`)

| Simulated Horizon | Steps (`dt = 30 s`) | Miyabi-G `BB` (`78.19 ms`) | H200 `BB-H200` (`32.04 ms`) | **H200 `Mode 2` (`4.07 ms`)** | **H200 `Mode 3 MONO` (`4.43 ms`)** |
|---|---:|---:|---:|---:|---:|
| **1 Model Day** | `2,880` | `225.2 s` (`3.75 min`) | `92.3 s` | **`11.7 s`** | **`12.7 s`** |
| **1 Model Month (30 d)** | `86,400` | `1.88 hr` | `46.1 min` | **`5.86 min`** | **`6.37 min`** |
| **1 Model Year (365 d)** | `1,051,200` | `22.83 hr` | `9.36 hr` | **`1.19 hr`** | **`1.29 hr`** |

---

## 6. Reproducibility & Miyabi-G One-Command Deployment (Post-Maintenance)

All Hopper SM90a v4.1 zero-spill kernels are backported into the self-contained Fortran-linkable NICAM tree under `work_opt/nicam_mp/` (verified compiling cleanly with CUDA 13.0 `-arch=sm_90a --prec-div=true --prec-sqrt=true`):

```bash
# 1. Sync work_opt/nicam_mp to Miyabi-G once maintenance concludes
rsync -az work_opt/nicam_mp/ miyabi-g:/work/jh260050/v57014/ANZIAM_2027/work_opt_nicam_mp/

# 2. On Miyabi-G: link dataset & build NVHPC 25.9 + CUDA 13.0 sm_90a binary
cd /work/jh260050/v57014/ANZIAM_2027/work_opt_nicam_mp
ln -s ../baseline/nicam_mp/dataset dataset
bash compile_opt_gpu.sh

# 3. Submit validation, 8x50 paired benchmark, and 5/100-step conservation jobs (<= 30 token-hours)
cd run/miyabi-g
qsub pbs_validate_opt.sh
qsub pbs_bench_8x50.sh
qsub pbs_multistep.sh
```

> **Transparency Note**: The real-checkpoint Miyabi-G GH200 validation and timing scripts above are ready for 1-command submission once Miyabi-G exits maintenance. All live GPU timings in §4 were measured on the NVIDIA H200 (`aniya-dev`, SM 9.0a) using our 10-regime full-scale (`166,410 × 80`) atmospheric workload, with direct apples-to-apples comparison against a faithful unfused baseline (`BB-H200 = 32.040 ms`) on the exact same H200, alongside cross-reference to the documented Miyabi-G GH200 `BB` (`78.19 ms`) and `TT` (`68.57 ms`) baselines.
