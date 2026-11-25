# NICAM NSW6 Cloud Microphysics: Structure-Preserving Mathematical Reformulation & Hopper GPU Acceleration

<p align="center">
  <strong><a href="#english">English</a></strong> &nbsp;|&nbsp; <strong><a href="#日本語">日本語</a></strong>
</p>

<p align="center">
  <img alt="Framework: Native CUDA C++17" src="https://img.shields.io/badge/Framework-Native%20CUDA%20C%2B%2B17-76B900?style=flat-square" />
  <img alt="Precision: Strict IEEE FP64" src="https://img.shields.io/badge/Precision-Strict%20IEEE%20FP64-00599C?style=flat-square" />
  <img alt="Architecture: NVIDIA Hopper SM90a" src="https://img.shields.io/badge/Target-Hopper%20SM90a%20%2F%20CUDA%20Graphs-76B900?style=flat-square" />
  <img alt="Register Spill: 0 Bytes" src="https://img.shields.io/badge/Register%20Spill-0%20Bytes%20(10%2F10%20Kernels)-2EA44F?style=flat-square" />
  <img alt="Validation: N_FAIL = 0" src="https://img.shields.io/badge/Fortran%20Oracle-N__FAIL%20%3D%200-2EA44F?style=flat-square" />
  <img alt="Speedup: 7.86x Same-GPU" src="https://img.shields.io/badge/H200%20Speedup-7.86%C3%97%20Same--GPU-E36209?style=flat-square" />
</p>

---

<a id="english"></a>
## English

### 1. Overview

This repository provides a **structure-preserving mathematical reformulation** and **native CUDA C++17 (`sm_90a`, CUDA Graphs) strict IEEE FP64 implementation** of the six-category bulk cloud microphysics scheme (**NSW6**) from the Non-hydrostatic Icosahedral Atmospheric Model (**NICAM**).

> **Implementation Framework**: This optimized implementation is written in **native NVIDIA CUDA C++17** (no OpenACC, Triton, or Kokkos runtime overhead), while isolating all mathematical reductions inside a portable, header-only `__host__ __device__` C++ core (`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`) so the same numerical formulations can also be evaluated inside C++/Kokkos NICAM.

Complementing algorithm-preserving Fortran-to-Kokkos migration efforts (Ren, Shimokawabe, & Yashiro, [`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313)), this work investigates how **algorithmic and mathematical restructuring** of the governing discrete microphysics—combined with spill-free GPU operator engineering—simultaneously achieves:

1. **Multiplicative runtime speedup** (`7.86×` on the same H200 GPU: `2.73×` GPU engineering × `2.88×` mathematical reduction), operating at `85%` of measured HBM3e streaming bandwidth;
2. **Strict IEEE double-precision fidelity** (`N_FAIL = 0` across all 23 output arrays / 34 checked diagnostics against the unmodified NICAM Fortran modules);
3. **Machine-epsilon multi-step conservation**, eliminating a severe trajectory-crossing defect in semi-Lagrangian sedimentation (`52.7% → 1.83 × 10⁻¹⁵` relative water-mass drift over 100 steps).

---

### 2. Research Goal, Cross-Framework Comparison & Todo

- **Goal**: Following **Ren et al. ([`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313))** as a benchmarking and evaluation guideline, perform a systematic, apples-to-apples performance and accuracy comparison across GPU programming models for NICAM NSW6—comparing the **Fortran/OpenACC** baseline, the **Triton** version, the **C++/Kokkos** version, and this repository's **Native CUDA C++** implementation (both algorithm-preserving and mathematically reformulated modes).
- **Todo (`@QIAN JIAWEI`)**:
  - [ ] **Check the performance of the Kokkos version** (`arXiv:2610.41313`) under the same grid size, real checkpoint, and timing protocol.
  - [ ] **Execute unified Miyabi-G (GH200) cross-framework benchmarks** comparing **OpenACC**, **Triton**, **Kokkos**, and **Native CUDA C++** (`Mode 1 ENG` vs. `Mode 2 ENG+MATH` vs. `Mode 3 MONO`) using Ren et al.'s evaluation methodology once Miyabi-G maintenance concludes.
  - [ ] **Evaluate portability of `nsw6_math_core.cuh` inside Kokkos** to cleanly decouple *programming-model overhead* (OpenACC / Triton / Kokkos / CUDA) from *mathematical-algorithmic speedup*.

#### Cross-Framework Comparison Matrix (Current Status)

| Implementation / Backend | Programming Model | Algorithm Variant | Measured Runtime (Full `166,410×80`) | Status |
| :--- | :--- | :--- | :---: | :--- |
| **NICAM Baseline** | Fortran + OpenACC | Original NSW6 | `78.19 ms` (Miyabi-G GH200) | Recorded baseline |
| **Triton Port** | Python / Triton | Original NSW6 | `68.57 ms` (Miyabi-G GH200) | Recorded baseline |
| **Kokkos Port** (`arXiv:2610.41313`) | C++ / Kokkos | Original NSW6 | **TBD** *(Miyabi-G GH200)* | **Todo (`@QIAN JIAWEI`)** |
| **This Repo — Mode 0 (`BB-H200`)** | Native CUDA C++ (unfused) | Original NSW6 | `31.90 ms` (NVIDIA H200) | Verified (`N_FAIL = 0`) |
| **This Repo — Mode 1 (`ENG`)** | Native CUDA C++ (Gen-4) | Original NSW6 | `11.70 ms` (NVIDIA H200) | Verified (`2.73×` vs. Mode 0) |
| **This Repo — Mode 2 (`ENG+MATH`)** | Native CUDA C++ + Graphs | Math-reformulated (exact) | **`4.06 ms`** / `3.89 ms` w/o `rceff` (H200) | Verified (`7.86×` / `8.21×`; GH200 PBS ready) |
| **This Repo — Mode 3 (`MONO`)** | Native CUDA C++ + Graphs | + Monotone sedimentation | **`4.41 ms`** (NVIDIA H200) | Verified (`7.23×`, `1.83×10⁻¹⁵` mass drift) |

> *Note*: Direct like-for-like cross-framework ratios must be measured on identical hardware (Miyabi-G GH200) and identical real-atmosphere checkpoints; H200 and GH200 numbers above are listed separately for transparency.

---

### 3. Key Results at a Glance

All GPU timings are measured on a single **NVIDIA H200 (141 GB HBM3e, CUDA 13.0)** at full NICAM production array size (**10 regions × 16,641 columns × 80 levels = 13.31M cells**, `dt = 30 s`), strictly compiled with `-O3 -arch=sm_90a --prec-div=true --prec-sqrt=true` (no `--use_fast_math`), using **8 cyclically rotated paired blocks** (95% Student-*t* confidence intervals, df = 7).

#### 3.1 Same-GPU Performance Attribution (H200, Full 13.31M Cells)

| Configuration | Description | Time / Step (ms) | Speedup vs. Mode 0 | Incremental Gain |
| :--- | :--- | :---: | :---: | :---: |
| **Mode 0 (`BB-H200`)** | Unfused CUDA transcription of original NSW6 formulas | `31.8989 ± 0.0058` | `1.00×` | — |
| **Mode 1 (`ENG`)** | Gen-4 Hopper kernels with original formulas | `11.7027 ± 0.0027` | `2.73×` | **`2.726×`** (GPU engineering) |
| **Mode 2 (`ENG+MATH`)** | + Mathematical reductions (exact, CUDA Graph, w/ `rceff`) | **`4.0575 ± 0.0027`** | **`7.86×`** | **`2.884×`** (Math reductions) |
| **Mode 2 (no `rceff`)** | Exact mode without optional effective-radius diagnostic | **`3.8858 ± 0.0023`** | **`8.21×`** | — |
| **Mode 3 (`ENG+MATH+MONO`)** | + $O(K)$ monotone backward-envelope sedimentation | **`4.4134 ± 0.0017`** | **`7.23×`** | `+0.356 ms` (`+8.8%`) |

- **Roofline efficiency**: Mode 2 sustains **`3.4 TB/s`** effective memory throughput (`85%` of the measured `4.02 TB/s` streaming ceiling; roofline lower bound `3.38 ms`).
- **Model-year impact**: Over one model year (`1,051,200` steps), microphysics runtime drops from **`9.31 h` to `1.19 h`**.
- **Cross-hardware orientation** *(non-like-for-like; different GPU & inputs)*: Compared against historical GH200 measurements of original Fortran/OpenACC (`78.19 ms`) and Triton (`68.57 ms`), Mode 2 corresponds to `19.3×` and `16.9×` ratios, respectively. Ready-to-run Miyabi-G GH200 PBS scripts are included under `work_opt/nicam_mp/run/miyabi-g/`.

#### 3.2 Work-Counter Reductions (Mode 1 → Mode 2, One Full Step)

| Operation / Metric | Reduction Factor | Source |
| :--- | :---: | :--- |
| `pow()` evaluations | **`131.5×` fewer** | Dyadic-lattice snow moments |
| `exp()` evaluations | **`7.31×` fewer** | Cancellation-free `expm1`-doubling |
| `log()` evaluations | **`5.73×` fewer** | Shared logarithm + dyadic recurrence |
| Saturation-adjustment iterations | **`2.46×` fewer** | Associated-Laguerre Householder-4 |
| Interior cells in heavy microphysics | **`6.0×` fewer** | Clear-sky dry-cell screening (`83.3%` skipped) |

#### 3.3 Validation & Multi-Step Conservation (`dt = 30 s`)

| Experiment | Modes | Water-Mass Residual | Moist Energy Residual | Notes |
| :--- | :--- | :---: | :---: | :--- |
| **Single-step Fortran oracle** (`1,024 × 80`) | Mode 1, 2 | — | — | **23 arrays / 34 checks: `N_FAIL = 0`**; `max\|ΔT\| = 3.2×10⁻¹⁰ K`, zero negatives |
| **Full-size, 5 steps (`150 s`)** | Mode 1, 2 | `6.2 × 10⁻²` | `1.9 × 10⁻⁴` | Non-monotone clipping breaks conservation |
| **Full-size, 5 steps (`150 s`)** | **Mode 3** | **`3.3 × 10⁻¹⁵`** | **`3.1 × 10⁻¹¹`** | `133k` full / `174k` half crossings removed |
| **Full-size, 100 steps (`3,000 s`)** | Mode 1, 2 | `5.3 × 10⁻¹` (`52.7%`) | `2.3 × 10⁻⁴` | Surface precipitation artificially inflated **`4.7×`** |
| **Full-size, 100 steps (`3,000 s`)** | **Mode 3** | **`1.8 × 10⁻¹⁵`** | **`2.4 × 10⁻¹⁰`** | **IEEE FP64 machine $\varepsilon$ (`14` orders better)**; `694k / 405k` crossings removed |
| **Real Fortran oracle, 100 steps** | Mode 0, 1, 2 | `8.9 × 10⁻²` | `1.4 × 10⁻³` | Linked directly against 20 NICAM Fortran modules |
| **Real Fortran oracle, 100 steps** | **Mode 3** | **`1.0 × 10⁻¹³`** | **`7.6 × 10⁻⁷`** | **`1.8 × 10³×`** energy residual reduction |

---

### 4. What Was Done

#### 4.1 Structure-Preserving Mathematical Reformulations

1. **Multi-Species Semi-Lagrangian Sedimentation (Fubini Factorisation & Monotone Projection)**
   - *Fubini summation reordering*: Reversing the species–level accumulation order collapses **28 species–attribute flux sweeps** into a single shared mass weight per column interval, eliminating redundant lookups and memory traffic.
   - *Exact monotonicity criterion & $O(K)$ backward-envelope projection*: An exact integral identity establishes that positivity and $L^1$-contraction hold **if and only if** interface departure points are non-decreasing. Under sharp precipitation tops at high Courant numbers, original trajectories cross, creating negative mixing ratios whose subsequent clipping breaks conservation. A minimal upper-constrained $O(K)$ backward-envelope projection removes all crossings while preserving mass to machine precision.
2. **Roh–Satoh (2014) Snow PSD Moments on the Dyadic Lattice $\frac{1}{8}\mathbb{Z}$**
   - All seven required snow moment orders lie on the dyadic lattice $\frac{1}{8}\mathbb{Z}$. Combining a single logarithm with cancellation-free `expm1`-doubling collapses **16 transcendental evaluations (`8 exp + 8 pow`) per cell** into **1 `log` + 2 `exp`** plus elementary products/square roots (maximum relative error `2.68 × 10⁻¹⁴`).
3. **Associated-Laguerre Structure of Saturation Adjustment**
   - Exact derivatives of the saturation mixing ratio satisfy the closed-form associated-Laguerre identity $q_{\mathrm{sat}}^{(m)}(T) = (-1)^m m!\, T^{-m} L_m^{(-\gamma)}(\alpha/T)\, q_{\mathrm{sat}}(T)$. This repairs a `~5%` Jacobian omission in NICAM's legacy Newton solver and yields a **quartic Householder-4 iteration** converging to double precision in at most 4 iterations.
4. **Clear-Sky Screening & Column Sparsity**
   - Over `83.3%` of interior atmospheric cells are clear and dry, and only `40.0%` of columns contain precipitation. Decoupling dry cells and non-precipitating columns removes severe warp divergence from heavy microphysics and sedimentation kernels.

#### 4.2 Four-Generation Hopper (`sm_90a`) CUDA Engineering

- **Gen 1 → Gen 2 (`9.12 ms → 5.71 ms`)**: Split sedimentation into half-level and full-level passes (eliminating `140 KB/block` shared-memory occupancy collapse); relocated physics constants to `__constant__` memory; eliminated `960 MB` of dead stores; fused saturation adjustment, negative filter, diagnostics, and cloud fraction into a single kernel.
- **Gen 3 (`5.71 ms → 5.34 ms`)**: Introduced a 44-register screening kernel streaming at `3.64 TB/s` that completes dry cells immediately and compacts the `16.7%` active wet cells via **warp-aggregated atomics** with packed 32-bit coordinates, enabling full-warp execution inside the 194-register microphysics kernel.
- **Gen 4 (`5.34 ms → 4.06 ms`)**: Relocated the 32-pointer state descriptor into `__constant__` memory (freeing `~64` registers/thread), achieving **0 bytes stack frame and 0 bytes register spill across all 10 kernels** (down from up to `752 B/thread` in Gen 1). Replaced persistent work queues with upper-bound early-exit grids for **full-step CUDA Graph capture**, added coalesced 3D vertical kinetic energy initialization, and skipped zero-density divisions and unchanged dry-level stores.

#### 4.3 Software Architecture & Portability

Rather than relying on high-level DSLs (OpenACC/Triton/PyTorch), the codebase uses a lightweight, zero-dependency three-layer architecture:

- **Native CUDA C++17 (`sm_90a` + CUDA Graphs)**: Full low-level control over register allocation, `__constant__` descriptors, warp-level stream compaction, and strict IEEE FP64 semantics (`work_opt/hopper_cuda/nsw6_hopper_v4.cuh`).
- **Portable Header-Only Math Core (`__host__ __device__`)**: All three mathematical reformulations live in a standalone C++ header (`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`), enabling identical CPU/GPU compilation and straightforward integration into C++/Kokkos NICAM codebases (`arXiv:2610.41313`).
- **Fortran 2003 `ISO_C_BINDING` Bridge**: Links directly against unmodified NICAM Fortran modules (`work_opt/nicam_mp/cuda_kernels/nsw6_cuda_kernels.cu`, `nsw6_cpu_bridge.cpp`).

---

### 5. Repository Layout

```text
.
├── baseline/                         # Pristine NICAM NSW6 Fortran/OpenACC baseline (UNMODIFIED, SHA-256 verified)
│   ├── SHA256SUMS                    # Authoritative checksums for baseline/nicam_mp/
│   └── nicam_mp/                     # Original NICAM microphysics mini-app
├── work_opt/
│   ├── hopper_cuda/                  # Standalone Hopper SM90a CUDA suite (v1–v4.1) & paired CI benchmarks
│   │   ├── nsw6_hopper_v4.cuh        # Final Gen-4 spill-free CUDA operators (CUDA Graph compatible)
│   │   ├── bench_paper_ci.cu         # 8-block paired timing + 95% Student-t CI + physical diagnostics
│   │   └── bench_h200_v3.cu          # Multi-step (5 / 100 steps) conservation driver
│   ├── nicam_mp/                     # Fortran-linkable optimized NICAM-MP tree
│   │   ├── cuda_kernels/
│   │   │   ├── nsw6_math_core.cuh    # Portable header-only CPU/GPU math core (3 theorems)
│   │   │   ├── nsw6_cuda_kernels.cu  # CUDA C++ entry point linked to Fortran via ISO_C_BINDING
│   │   │   └── nsw6_cpu_bridge.cpp   # Host C++ bridge for Fortran differential oracle
│   │   ├── compile_opt_gpu.sh        # Miyabi-G GH200 NVHPC + CUDA 13.0 build script
│   │   └── run/miyabi-g/             # Ready-to-run PBS validation, 8×50 benchmark & multi-step scripts
│   └── analysis/
│       ├── verify_fortran_oracle.f90 # Differential oracle linked against 20 real NICAM Fortran modules
│       └── verify_math_and_kernels.cpp
├── docs/
│   └── ANZIAM_2027_REPORT.md         # Detailed technical engineering report
└── output/                           # Reference project brief & documents
```

---

### 6. Build, Verification & Reproducibility

#### 6.1 Verify Baseline Integrity
```bash
cd baseline && shasum -a 256 -c SHA256SUMS
```

#### 6.2 Run Full-Scale Hopper H200 Paired CI Benchmark (Strict IEEE FP64)
```bash
cd work_opt/hopper_cuda
nvcc -O3 -std=c++17 -arch=sm_90a --prec-div=true --prec-sqrt=true \
     -Xptxas -v bench_paper_ci.cu -o bench_paper_ci
./bench_paper_ci
```

#### 6.3 Build & Run on Miyabi-G (GH200)
```bash
cd work_opt/nicam_mp
./compile_opt_gpu.sh
qsub run/miyabi-g/pbs_validate_opt.sh   # Real-checkpoint differential validation
qsub run/miyabi-g/pbs_bench_8x50.sh     # 8-block paired benchmark
qsub run/miyabi-g/pbs_multistep.sh      # Multi-step budget verification
```

---

<a id="日本語"></a>
## 日本語

### 1. 概要

本リポジトリは、全球非静力学正二十面体大気モデル **NICAM** の 6 カテゴリー・バルク雲微物理スキーム（**NSW6**）に対して、**構造保存型の数理再定式化**および **ネイティブ CUDA C++17（Hopper `sm_90a` / CUDA グラフ）による厳密倍精度（IEEE FP64）実装**を提供するものです。

> **実装フレームワーク**：本リポジトリの最適化版は **ネイティブ NVIDIA CUDA C++17** で実装されています（OpenACC / Triton / Kokkos 等のランタイムを介さずハードウェアを直接制御）。一方、3 つの数理再定式化は外部依存のない Header-only `__host__ __device__` C++ コア（`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`）として分離しており、C++/Kokkos 版 NICAM への移植・比較評価も容易に行える設計です。

既存アルゴリズムを忠実に保持した Fortran→Kokkos 移植研究（Ren, Shimokawabe, & Yashiro, [`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313)）の次段階として、離散物理スキーム内部の数理構造そのものを再構成し、スピルゼロの GPU オペレータ設計と組み合わせることで、以下の 3 点を同時に達成しました。

1. **実装最適化と数理改革の乗算的な高速化**（同一 H200 GPU 上で **`7.86 倍`**：GPU 実装最適化 `2.73 倍` × 数理再定式化 `2.88 倍`、実効メモリ帯域幅は実測ストリーム帯域の `85%` に到達）
2. **厳密倍精度（IEEE FP64）での数値的整合性**（元 NICAM Fortran モジュール 20 本との差分検証において全 23 出力配列・34 項目で **`N_FAIL = 0`**）
3. **多数ステップ積算におけるマシンイプシロン級の保存性回復**（半ラグランジュ鉛直沈降における軌跡交差欠陥を解消し、100 ステップ連続積分時の水質量相対誤差を **`52.7% → 1.83 × 10⁻¹⁵`** へ劇的に低減）

---

### 2. 研究目標・フレームワーク横断比較・今後の課題（Todo）

- **目標（Goal）**：**Ren らの論文（[`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313)）の評価指針（ガイドライン）** に従い、NICAM NSW6 の各 GPU 実装バックエンド——**Fortran/OpenACC 版**、**Triton 版**、**C++/Kokkos 版**、および本リポジトリの **ネイティブ CUDA C++ 版**（アルゴリズム保存モードおよび数理再定式化モード）——の性能と数値精度を同一条件で体系的に比較評価する。
- **Todo（`@QIAN JIAWEI`）**：
  - [ ] **Kokkos 版（`arXiv:2610.41313`）の性能計測・検証**：同一格子規模・実チェックポイント・統一計測プロトコルにて Kokkos 版の実行時間を測定する。
  - [ ] **Miyabi-G（GH200）での横断ベンチマーク実施**：メンテナンス完了後、Ren らのガイドラインに基づき **OpenACC / Triton / Kokkos / ネイティブ CUDA C++（Mode 1 / 2 / 3）** の同一ハードウェア比較を完了させる。
  - [ ] **Kokkos バックエンドへの `nsw6_math_core.cuh` 統合検証**：フレームワーク起因のオーバーヘッド（OpenACC vs. Triton vs. Kokkos vs. CUDA）と数理アルゴリズム改革による高速化を明確に分離する。

#### バックエンド別比較マトリクス（現在の状況）

| 実装 / バックエンド | フレームワーク | アルゴリズム | 実測実行時間（全領域 `166,410×80`） | ステータス |
| :--- | :--- | :--- | :---: | :--- |
| **NICAM 元ベースライン** | Fortran + OpenACC | 原 NSW6 | `78.19 ms`（Miyabi-G GH200） | 計測済ベースライン |
| **Triton 移植版** | Python / Triton | 原 NSW6 | `68.57 ms`（Miyabi-G GH200） | 計測済ベースライン |
| **Kokkos 移植版**（`arXiv:2610.41313`） | C++ / Kokkos | 原 NSW6 | **未計測（TBD）** *(Miyabi-G GH200)* | **Todo（`@QIAN JIAWEI`）** |
| **本リポジトリ — Mode 0 (`BB-H200`)** | ネイティブ CUDA C++（非融合） | 原 NSW6 | `31.90 ms`（NVIDIA H200） | 検証済（`N_FAIL = 0`） |
| **本リポジトリ — Mode 1 (`ENG`)** | ネイティブ CUDA C++（第 4 世代） | 原 NSW6 | `11.70 ms`（NVIDIA H200） | 検証済（Mode 0 比 `2.73×`） |
| **本リポジトリ — Mode 2 (`ENG+MATH`)** | ネイティブ CUDA C++ + グラフ | 数理再定式化（厳密等価） | **`4.06 ms`** / `rceff` なし `3.89 ms`（H200） | 検証済（`7.86×` / `8.21×`、GH200 用 PBS 準備済） |
| **本リポジトリ — Mode 3 (`MONO`)** | ネイティブ CUDA C++ + グラフ | ＋単調後向包絡射影 | **`4.41 ms`**（NVIDIA H200） | 検証済（`7.23×`、質量誤差 `1.83×10⁻¹⁵`） |

> *注*：厳密なフレームワーク間比較は同一ハードウェア（Miyabi-G GH200）および同一実データ入力で行う必要があるため、上記表では H200 と GH200 の数値を明確に区別して記載しています。

---

### 3. 主な成果サマリー

GPU 性能計測は **NVIDIA H200（141 GB HBM3e, CUDA 13.0）** 1 基上にて、NICAM 本番相当の全領域格子（**10 リージョン × 16,641 鉛直柱 × 80 層 = 1,331 万格子点**、`dt = 30 s`）を対象に実施しました。全カーネルを厳密倍精度オプション（`-O3 -arch=sm_90a --prec-div=true --prec-sqrt=true`、`--use_fast_math` 禁止）でコンパイルし、**実行順序を巡回させた 8 ペアブロック計測**（自由度 7 の Student-*t* 95% 信頼区間）により統計評価しています。

#### 3.1 同一 GPU における要因分解（NVIDIA H200、全 1,331 万格子点）

| 構成 | 内容 | 1 ステップ実行時間 (ms) | 対 Mode 0 速度比 | 段階別寄与 |
| :--- | :--- | :---: | :---: | :---: |
| **Mode 0 (`BB-H200`)** | 原 NSW6 数式をそのまま展開した非融合 CUDA ベースライン | `31.8989 ± 0.0058` | `1.00×` | — |
| **Mode 1 (`ENG`)** | 第 4 世代 Hopper カーネル＋原数式（実装最適化のみ） | `11.7027 ± 0.0027` | `2.73×` | **`2.726×`**（GPU 実装） |
| **Mode 2 (`ENG+MATH`)** | ＋数理再定式化（厳密等価、CUDA グラフ、`rceff` 含む） | **`4.0575 ± 0.0027`** | **`7.86×`** | **`2.884×`**（数理改革） |
| **Mode 2 (`rceff` なし)** | 任意の有効半径診断を除外した厳密等価モード | **`3.8858 ± 0.0023`** | **`8.21×`** | — |
| **Mode 3 (`ENG+MATH+MONO`)** | ＋$O(K)$ 単調後向包絡射影（保存性修復モード） | **`4.4134 ± 0.0017`** | **`7.23×`** | `+0.356 ms`（`+8.8%`） |

- **ルーフライン到達度**：Mode 2 の実効メモリ帯域幅は **`3.4 TB/s`** に達し、同一 GPU の単純ストリームカーネル実測値（`4.02 TB/s`）の **`85%`**（理論下限 `3.38 ms` に肉薄）を達成しています。
- **気候シミュレーション換算**：1 モデル年（`1,051,200` ステップ）あたりの雲微物理計算時間は **`9.31 時間` → `1.19 時間`** に短縮されます。
- **参考値（異種 GPU 比較）**：Miyabi-G（GH200）上の既存 Fortran/OpenACC 版（`78.19 ms`）および Triton 版（`68.57 ms`）と比較すると、Mode 2 はそれぞれ `19.3 倍`・`16.9 倍` に相当します（※ハードウェア・入力条件が異なるため参考値）。Miyabi-G 向けの実データ検証・ベンチマーク用 PBS スクリプトは `work_opt/nicam_mp/run/miyabi-g/` に完備しています。

#### 3.2 演算量・処理格子数の削減（Mode 1 → Mode 2、1 ステップ全体）

| 指標 | 削減倍率 | 主な要因 |
| :--- | :---: | :--- |
| `pow()` 呼び出し回数 | **`131.5 分の 1`** | 二進格子上の雪粒子モーメント評価 |
| `exp()` 呼び出し回数 | **`7.31 分の 1`** | 桁落ち回避型 `expm1` 倍増漸化式 |
| `log()` 呼び出し回数 | **`5.73 分の 1`** | 対数の共通化と二進漸化計算 |
| 飽和調節の反復回数 | **`2.46 分の 1`** | 陪ラゲール多項式による 4 次ハウスホルダー法 |
| 重微物理カーネルの処理格子点 | **`6.0 分の 1`** | 晴天・乾燥格子の事前スクリーニング（`83.3%` を除外） |

#### 3.3 差分検証および多数ステップ保存性（`dt = 30 s`）

| 検証実験 | モード | 水質量相対残差 | 湿潤全エネルギー相対残差 | 備考 |
| :--- | :--- | :---: | :---: | :--- |
| **単ステップ Fortran 差分検証**（`1,024 × 80`） | Mode 1, 2 | — | — | **全 23 配列・34 項目で `N_FAIL = 0`**、最大 $\|\Delta T\| = 3.2\times 10^{-10}\,\mathrm{K}$、負値 0 |
| **全領域・5 ステップ（`150 s`）** | Mode 1, 2 | `6.2 × 10⁻²` | `1.9 × 10⁻⁴` | 非単調クリッピングによる保存性破綻 |
| **全領域・5 ステップ（`150 s`）** | **Mode 3** | **`3.3 × 10⁻¹⁵`** | **`3.1 × 10⁻¹¹`** | 軌跡交差 `133k`（全層）/ `174k`（半層）を除去 |
| **全領域・100 ステップ（`3,000 s`）** | Mode 1, 2 | `5.3 × 10⁻¹`（`52.7%`） | `2.3 × 10⁻⁴` | 累積地表面降水量が **`4.7 倍`** に過大評価 |
| **全領域・100 ステップ（`3,000 s`）** | **Mode 3** | **`1.8 × 10⁻¹⁵`** | **`2.4 × 10⁻¹⁰`** | **FP64 機械精度（14 桁改善）**、交差 `694k / 405k` を除去 |
| **実 Fortran オラクル・100 ステップ** | Mode 0, 1, 2 | `8.9 × 10⁻²` | `1.4 × 10⁻³` | 原 NICAM Fortran 20 モジュールと直接リンク |
| **実 Fortran オラクル・100 ステップ** | **Mode 3** | **`1.0 × 10⁻¹³`** | **`7.6 × 10⁻⁷`** | エネルギー残差を約 **`1.8 × 10³ 倍`** 改善 |

---

### 4. 実施した最適化内容

#### 4.1 構造保存型の数理再定式化

1. **多水物質・半ラグランジュ鉛直沈降（フビニ因数分解と単調後向包絡射影）**
   - *フビニ和順序交換*：水物質種・物理量ごとの **28 系統の鉛直フラックス走査**を「層優先・種後置」へ組み替え、単一の共有質量重みに集約することで重複メモリアクセスを大幅に削減しました。
   - *厳密単調性条件と $O(K)$ 後向包絡射影*：厳密積分恒等式により、非負性および $L^1$ 縮小性の必要十分条件が「界面出発点の単調非減少性」であることを証明しました。強降水頂部・高クーラン数条件下で生じる出発点交差と後続クリッピングによる質量不保存を、$O(K)$ の最小上界制約・後向包絡射影により完全に除去します。
2. **Roh–Satoh (2014) 雪粒子粒径分布モーメントの二進格子 $\frac{1}{8}\mathbb{Z}$ 構造**
   - 必要な 7 つのモーメント次数がすべて二進格子 $\frac{1}{8}\mathbb{Z}$ に属することに着目し、桁落ちを避ける `expm1` 倍増漸化式により、1 格子点あたり **16 回の超越関数（`8 exp + 8 pow`）** を **1 回の `log` と 2 回の `exp`** に圧縮しました（最大相対誤差 `2.68 × 10⁻¹⁴`）。
3. **飽和調節（Saturation Adjustment）の陪ラゲール構造と 4 次収束解法**
   - 飽和混合比の任意階温度微分が陪ラゲール多項式 $q_{\mathrm{sat}}^{(m)}(T) = (-1)^m m!\, T^{-m} L_m^{(-\gamma)}(\alpha/T)\, q_{\mathrm{sat}}(T)$ で厳密に表せることを導出し、既存 NICAM ニュートン法のヤコビアン欠落（約 `5%` の誤差）を修正。最大 4 反復で倍精度収束する **4 次ハウスホルダー法**へ高度化しました。
4. **乾燥格子スクリーニングと降水柱の疎性活用**
   - 大気内部の `83.3%` を占める晴天乾燥格子を軽量カーネルで即座に完了させ、`16.7%` の湿潤格子のみを重微物理カーネルへ送出。さらに沈降処理も降水を含む `40.0%` の鉛直柱のみに限定しました。

#### 4.2 Hopper（`sm_90a`）向け 4 世代 CUDA カーネル最適化

- **Gen 1 → Gen 2（`9.12 ms → 5.71 ms`）**：沈降を半層・全層の 2 パスに分離し、`140 KB/block` の共有メモリ占有によるオキュパンシ低下を解消。物理定数の `__constant__` 化、`960 MB` の冗長ストア除去、さらに飽和調節・負値フィルタ・診断量・雲量の 1 カーネル融合を実施。
- **Gen 3（`5.71 ms → 5.34 ms`）**：`3.64 TB/s` で動作する 44 レジスタのスクリーニングカーネルを導入し、32-bit 圧縮座標と **ワープ集約アトミック（warp-aggregated atomics）** で湿潤格子キューを構築。194 レジスタの重微物理カーネルをワープ分岐なしで高効率実行。
- **Gen 4（`5.34 ms → 4.06 ms`）**：32 本の状態ポインタ記述子も `__constant__` メモリへ移設し（スレッドあたり約 64 レジスタ削減）、**全 10 カーネルでスタックフレーム 0 B・レジスタスピル 0 B** を達成（Gen 1 では最大 `752 B/thread`）。永続キューを上限サイズ＋早期終了グリッドに置換して **1 ステップ全体を CUDA グラフ化**し、3 次元コアレス鉛直運動エネルギー初期化およびゼロ除算・乾燥層再書込の回避を実装。

#### 4.3 採用フレームワークとソフトウェア構成

高水準 DSL（OpenACC / Triton / PyTorch 等）には依存せず、厳密倍精度科学計算と移植性を両立する 3 層構成を採用しています。

- **ネイティブ CUDA C++17（`sm_90a` + CUDA グラフ）**：レジスタ割当、`__constant__` 記述子、ワープ集約ストリーム圧縮、および厳密 IEEE FP64 演算を完全に制御（`work_opt/hopper_cuda/nsw6_hopper_v4.cuh`）。
- **移植性の高い Header-only 数理コア（`__host__ __device__` C++）**：3 つの数理再定式化を外部依存ゼロの単一ヘッダ（`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`）に集約。CPU/GPU 同一ソースでの検証に加え、`arXiv:2610.41313` の **C++/Kokkos 版 NICAM** への組み込みも容易です。
- **Fortran 2003 `ISO_C_BINDING` ブリッジ**：無改変の元 NICAM Fortran モジュール群と直接リンクして動作・差分検証可能（`work_opt/nicam_mp/cuda_kernels/nsw6_cuda_kernels.cu`, `nsw6_cpu_bridge.cpp`）。

---

### 5. ビルド・検証手順

```bash
# 1. 元ベースラインの無改変確認（SHA-256 チェック）
cd baseline && shasum -a 256 -c SHA256SUMS

# 2. NVIDIA H200 上での厳密倍精度 8 ペアブロック統計ベンチマーク
cd work_opt/hopper_cuda
nvcc -O3 -std=c++17 -arch=sm_90a --prec-div=true --prec-sqrt=true \
     -Xptxas -v bench_paper_ci.cu -o bench_paper_ci
./bench_paper_ci

# 3. Miyabi-G (GH200) 向けビルドおよび PBS 検証ジョブ投入
cd work_opt/nicam_mp
./compile_opt_gpu.sh
qsub run/miyabi-g/pbs_validate_opt.sh   # 実チェックポイント差分検証
qsub run/miyabi-g/pbs_bench_8x50.sh     # 8 ブロック性能評価
qsub run/miyabi-g/pbs_multistep.sh      # 多数ステップ保存性検証
```
