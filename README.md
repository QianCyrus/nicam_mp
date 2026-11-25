# NICAM NSW6 Cloud Microphysics — Native CUDA C++ & Mathematical Reformulation

<p align="center">
  <strong><a href="#english">English</a></strong> &nbsp;|&nbsp; <strong><a href="#日本語">日本語</a></strong>
</p>

<p align="center">
  <img alt="Framework: Native CUDA C++17" src="https://img.shields.io/badge/Framework-Native%20CUDA%20C%2B%2B17-76B900?style=flat-square" />
  <img alt="Precision: Strict IEEE FP64" src="https://img.shields.io/badge/Precision-Strict%20IEEE%20FP64-00599C?style=flat-square" />
  <img alt="Spill: 0 Bytes" src="https://img.shields.io/badge/Register%20Spill-0%20Bytes-2EA44F?style=flat-square" />
  <img alt="Oracle: N_FAIL = 0" src="https://img.shields.io/badge/Fortran%20Oracle-N__FAIL%20%3D%200-2EA44F?style=flat-square" />
  <img alt="H200 Speedup: 7.86x" src="https://img.shields.io/badge/H200%20Speedup-7.86%C3%97-E36209?style=flat-square" />
</p>

---

<a id="english"></a>
## English

### Overview & Goal
- **Framework**: **Native NVIDIA CUDA C++17** (`sm_90a`, CUDA Graphs, strict IEEE FP64 without `--use_fast_math`), with mathematical reductions isolated in a portable header-only C++ core (`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`) reusable in Kokkos.
- **Goal**: Using **Ren et al. ([`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313))** as a guideline, systematically compare NICAM NSW6 performance across **Fortran/OpenACC**, **Triton**, **C++/Kokkos**, and **Native CUDA C++** (both algorithm-preserving and mathematically reformulated modes).
- **Todo (`@QIAN JIAWEI`)**:
  - [ ] Check and benchmark the **Kokkos version** (`arXiv:2610.41313`) under the same domain and timing protocol.
  - [ ] Run unified **Miyabi-G (GH200)** real-checkpoint comparisons across **OpenACC / Triton / Kokkos / Native CUDA C++** post-maintenance.

### What Was Done
1. **Sedimentation (Fubini + Monotone Projection)**: Reordered 28 species–attribute flux sweeps into 1 shared mass weight; proved positivity/$L^1$-contraction $\iff$ monotone departure points and added an $O(K)$ backward-envelope projection removing all trajectory crossings.
2. **Snow Moments (Dyadic Lattice $\frac{1}{8}\mathbb{Z}$)**: Reduced 16 transcendentals (`8 exp + 8 pow`) per cell to `1 log + 2 exp` via cancellation-free `expm1`-doubling (`131.5×` fewer `pow`, `7.31×` fewer `exp`, `5.73×` fewer `log`; max rel. err `2.68×10⁻¹⁴`).
3. **Saturation Adjustment (Laguerre Householder-4)**: Repaired a `~5%` Jacobian omission using associated-Laguerre derivatives $L_m^{(-\gamma)}(\alpha/T)$, converging in `≤ 4` steps (`2.46×` fewer iterations).
4. **Hopper `sm_90a` CUDA Engineering (Gen 1→4)**: 44-register dry-cell screening (`83.3%` dry cells skipped, `6.0×` fewer active microphysics cells) + warp-aggregated atomic compaction; `__constant__` state descriptors achieving **0 B register spill across all 10 kernels**; full-step **CUDA Graph** reaching **`3.4 TB/s`** (`85%` of measured H200 streaming bandwidth).

### Performance & Conservation Summary (`166,410 × 80` = 13.31M cells, `dt = 30 s`)

| Backend / Mode | Framework | Algorithm | Time / Step | Notes / Status |
| :--- | :--- | :--- | :---: | :--- |
| **Baseline** | Fortran + OpenACC | Original NSW6 | `78.19 ms` *(GH200)* | Recorded baseline |
| **Triton** | Python / Triton | Original NSW6 | `68.57 ms` *(GH200)* | Recorded baseline |
| **Kokkos** (`arXiv:2610.41313`) | C++ / Kokkos | Original NSW6 | **TBD** *(GH200)* | **Todo (`@QIAN JIAWEI`)** |
| **Mode 0 (`BB-H200`)** | Native CUDA (unfused) | Original NSW6 | `31.90 ms` *(H200)* | `1.00×` same-GPU reference |
| **Mode 1 (`ENG`)** | Native CUDA (Gen-4) | Original NSW6 | `11.70 ms` *(H200)* | **`2.73×`** (engineering only) |
| **Mode 2 (`ENG+MATH`)** | Native CUDA + Graph | Math-reformulated | **`4.06 ms`** *(H200)* | **`7.86×`** (`2.73×` eng × `2.88×` math); `3.89 ms` (`8.21×`) w/o `rceff` |
| **Mode 3 (`MONO`)** | Native CUDA + Graph | + Monotone proj. | **`4.41 ms`** *(H200)* | **`7.23×`** (`+0.36 ms`); fixes multi-step mass conservation |

- **Validation**: `N_FAIL = 0` across all 23 arrays / 34 checks vs. 20 unmodified NICAM Fortran modules (`max |ΔT| = 3.2×10⁻¹⁰ K`, zero negatives).
- **100-Step Conservation (`3,000 s`)**: Mode 3 cuts relative water-mass drift from **`52.7% → 1.83×10⁻¹⁵`** (eliminating `4.7×` inflated surface precipitation) and moist energy drift from **`2.3×10⁻⁴ → 2.4×10⁻¹⁰`**.
- **Model-Year Impact**: `1,051,200` steps drop from **`9.31 h → 1.19 h`** on H200 (GH200 PBS scripts ready under `work_opt/nicam_mp/run/miyabi-g/`).

### Quick Start
```bash
cd baseline && shasum -a 256 -c SHA256SUMS                         # Verify untouched baseline
cd ../work_opt/hopper_cuda && nvcc -O3 -std=c++17 -arch=sm_90a \
   --prec-div=true --prec-sqrt=true bench_paper_ci.cu -o bench && ./bench  # H200 8-block CI
cd ../nicam_mp && ./compile_opt_gpu.sh                             # Miyabi-G GH200 build
```

---

<a id="日本語"></a>
## 日本語

### 概要・研究目標
- **実装フレームワーク**：**ネイティブ NVIDIA CUDA C++17**（`sm_90a`・CUDA グラフ・厳密倍精度 IEEE FP64、`--use_fast_math` 禁止）。数理コアは Kokkos にも流用可能な Header-only C++（`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`）として分離しています。
- **目標（Goal）**：**Ren らの論文（[`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313)）をガイドライン**として、NICAM NSW6 における **Fortran/OpenACC 版**・**Triton 版**・**C++/Kokkos 版**・本リポジトリの **ネイティブ CUDA C++ 版**の性能と精度を体系的に比較評価します。
- **Todo（`@QIAN JIAWEI`）**：
  - [ ] **Kokkos 版（`arXiv:2610.41313`）の性能確認・ベンチマーク計測**を実施する。
  - [ ] Miyabi-G（GH200）復旧後、**OpenACC / Triton / Kokkos / ネイティブ CUDA C++** の同一実チェックポイント比較を完了させる。

### 主な最適化内容
1. **鉛直沈降（フビニ因数分解＋単調射影）**：28 系統のフラックス走査を 1 つの共有質量重みに集約。さらに $O(K)$ 単調後向包絡射影により出発点交差を完全除去し、非負性・保存性を保証。
2. **雪モーメント（二進格子 $\frac{1}{8}\mathbb{Z}$）**：桁落ち回避 `expm1` 倍増により 1 格子点 `16` 回の超越関数を `1 log + 2 exp` に圧縮（`pow 131.5×`・`exp 7.31×`・`log 5.73×` 削減、最大相対誤差 `2.68×10⁻¹⁴`）。
3. **飽和調節（陪ラゲール 4 次ハウスホルダー法）**：既存ニュートン法のヤコビアン欠落（約 `5%`）を $L_m^{(-\gamma)}(\alpha/T)$ により厳密修復し、最大 4 反復で倍精度収束（反復数 `2.46×` 削減）。
4. **Hopper `sm_90a` CUDA 最適化（Gen 1→4）**：乾燥格子 `83.3%` の事前スクリーニング（微物理処理格子 `6.0×` 削減）＋ワープ集約アトミック圧縮、`__constant__` 記述子化による**全 10 カーネル・レジスタスピル 0 B**、1 ステップ全体 **CUDA グラフ化**で実効帯域 **`3.4 TB/s`**（実測上限の `85%`）を達成。

### 性能・保存性サマリー（全領域 `166,410 × 80` = 1,331 万格子点、`dt = 30 s`）

| バックエンド / モード | フレームワーク | アルゴリズム | 実行時間 / step | 備考・ステータス |
| :--- | :--- | :--- | :---: | :--- |
| **元ベースライン** | Fortran + OpenACC | 原 NSW6 | `78.19 ms` *(GH200)* | 計測済ベースライン |
| **Triton 版** | Python / Triton | 原 NSW6 | `68.57 ms` *(GH200)* | 計測済ベースライン |
| **Kokkos 版**（`arXiv:2610.41313`） | C++ / Kokkos | 原 NSW6 | **未計測** *(GH200)* | **Todo（`@QIAN JIAWEI`）** |
| **Mode 0 (`BB-H200`)** | ネイティブ CUDA（非融合） | 原 NSW6 | `31.90 ms` *(H200)* | 同一 GPU 基準（`1.00×`） |
| **Mode 1 (`ENG`)** | ネイティブ CUDA（Gen-4） | 原 NSW6 | `11.70 ms` *(H200)* | **`2.73×`**（実装最適化のみ） |
| **Mode 2 (`ENG+MATH`)** | ネイティブ CUDA + グラフ | 数理再定式化 | **`4.06 ms`** *(H200)* | **`7.86×`**（実装 `2.73×` × 数理 `2.88×`）、`rceff` なし `3.89 ms`（`8.21×`） |
| **Mode 3 (`MONO`)** | ネイティブ CUDA + グラフ | ＋単調射影 | **`4.41 ms`** *(H200)* | **`7.23×`**（`+0.36 ms`）、多数ステップ保存性修復 |

- **差分検証**：元 NICAM Fortran 20 モジュールに対し全 23 配列・34 項目 **`N_FAIL = 0`**（最大 $\|\Delta T\| = 3.2\times 10^{-10}\,\mathrm{K}$、負値 0）。
- **100 ステップ保存性（`3,000 s`）**：Mode 3 により水質量相対誤差 **`52.7% → 1.83×10⁻¹⁵`**（地表面降水の `4.7×` 過大評価を解消）、湿潤全エネルギー誤差 **`2.3×10⁻⁴ → 2.4×10⁻¹⁰`** へ改善。
- **リポジトリ主要構成**：`baseline/`（無改変原本・SHA-256 検証済）、`work_opt/hopper_cuda/`（CUDA 実装・ベンチマーク）、`work_opt/nicam_mp/`（Fortran 接続版・Miyabi-G 用 PBS スクリプト）、`work_opt/analysis/`（Fortran 差分検証）、`docs/`（詳細レポート）。
