# NICAM NSW6 Cloud Microphysics — Mathematical Reformulation on Native CUDA & Portable Kokkos

<p align="center">
  <strong><a href="#english">English</a></strong> &nbsp;|&nbsp; <strong><a href="#日本語">日本語</a></strong>
</p>

<p align="center">
  <img alt="Frameworks: Native CUDA C++17 & Kokkos 4.7" src="https://img.shields.io/badge/Frameworks-Native%20CUDA%20C%2B%2B17%20%7C%20Kokkos%204.7-76B900?style=flat-square" />
  <img alt="Precision: Strict IEEE FP64" src="https://img.shields.io/badge/Precision-Strict%20IEEE%20FP64-00599C?style=flat-square" />
  <img alt="Spill: 0 Bytes" src="https://img.shields.io/badge/Register%20Spill-0%20Bytes-2EA44F?style=flat-square" />
  <img alt="Oracle: N_FAIL = 0" src="https://img.shields.io/badge/Validation-N__FAIL%20%3D%200-2EA44F?style=flat-square" />
  <img alt="H200 Speedup: 7.85x CUDA / 5.81x Kokkos" src="https://img.shields.io/badge/H200%20Speedup-7.85%C3%97%20CUDA%20%7C%205.81%C3%97%20Kokkos-E36209?style=flat-square" />
</p>

---

<a id="english"></a>
## English

### Overview, Goal & Todo
- **Frameworks**: Implemented in both **Native NVIDIA CUDA C++17** (`sm_90a`, CUDA Graphs) and **100% Pure Portable Kokkos 4.7** (`View` / `RangePolicy` / `TeamPolicy` scratch compaction, zero raw CUDA syntax), sharing one header-only C++ math core (`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`), compiled in strict IEEE FP64 (`--prec-div=true --prec-sqrt=true`, no `--use_fast_math`, **0 B register spill** across all kernels).
- **Goal**: Following **Ren et al. ([`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313))** as a guideline, compare **Fortran/OpenACC**, **Triton**, **C++/Kokkos**, and **Native CUDA C++** for NICAM NSW6.
- **Todo (`@QIAN JIAWEI`)**:
  - [x] **Benchmark Kokkos on idle NVIDIA H200**: profiled upstream [`Masahiro-Ren/nicam_kokkos_microphysics`](https://github.com/Masahiro-Ren/nicam_kokkos_microphysics) (`71.05 ms` wall-clock / `20.28 ms` pure GPU kernels for 1 region `16,641×96`) and ran full-scale (`166,410×80`) 8-block paired CI tests in **100% pure portable Kokkos 4.7** (**`28.71 → 13.64 → 4.94 ms`**, **`5.81×`** total / **`2.76×`** pure math, `N_FAIL = 0`).
  - [ ] Run unified **Miyabi-G (GH200)** real-checkpoint comparisons across **OpenACC / Triton / Kokkos / Native CUDA C++** post-maintenance.

### Mathematical Methods & Core Gains
1. **Snow Moments (Dyadic Lattice $n \in \frac{1}{8}\mathbb{Z}$) + Dry Screening**: Cancellation-free `expm1`-doubling reduces 16 transcendentals (`8 exp + 8 pow`) per cell to `1 log + 2 exp` (`131.5×` fewer `pow`, `7.31×` fewer `exp`, `5.73×` fewer `log`), plus `83.3%` clear-sky dry-cell screening (`6.0×` fewer active microphysics cells).
2. **Sedimentation (Fubini Factorisation + $O(K)$ Monotone Projection)**: Shared mass weight $W_{k,j}=\sum_s L_s(k,j)a_{j,s}$ collapses 28 flux sweeps on `40%` precipitating columns; backward-envelope projection $\tilde{\xi}_k=\min(\xi_k,\tilde{\xi}_{k+1})$ (`+0.35–0.36 ms`) eliminates all trajectory crossings and cuts 100-step water drift from **`52.7%` → machine $\varepsilon$**.
3. **Saturation Adjustment (Associated-Laguerre Householder-4)**: Exact derivatives $q_{\mathrm{sat}}^{(m)}(T)=(-1)^m m!\,T^{-m} L_m^{(-\gamma)}(\alpha/T)\,q_{\mathrm{sat}}(T)$ fix a `~5%` Jacobian defect, converging in `≤ 4` steps (`2.46×` fewer iterations).

### Idle-H200 Same-Session Comparison (`166,410 × 80` = 13.31M cells, `dt = 30 s`, 8 Paired Blocks, 95% CI)

| Implementation / Mode | Framework | Algorithm | Step Time (ms) | Speedup & Validation |
| :--- | :--- | :--- | :---: | :--- |
| **Baseline** | Fortran + OpenACC | Original NSW6 | `78.19` *(GH200)* | Historical GH200 baseline |
| **Triton** | Python / Triton | Original NSW6 | `68.57` *(GH200)* | Historical GH200 baseline |
| **Ren Kokkos (upstream)** | Kokkos 4.7 | Original NSW6 (`16,641×96`, 1 reg) | `71.05` wall / `20.28` GPU *(H200)* | 889 launches + 228 alloc/free (`~203 ms` 10-reg GPU eq.) |
| **Kokkos-BB** | Pure Kokkos 4.7 | Original NSW6 (unfused, 0 alloc) | `28.71 ± 0.09` *(H200)* | `1.00×` Kokkos reference |
| **Kokkos-ENG** | Pure Kokkos 4.7 | Original NSW6 (whole-domain) | `13.64 ± 0.14` *(H200)* | **`2.11×`** (Kokkos eng); **`N_FAIL = 0`** |
| **Kokkos-MATH** *(w/ `rceff`)* | **Pure Kokkos 4.7** | **+ Our 3 Math Theorems** | **`4.94 ± 0.02`** *(H200)* | **`5.81×`** (`2.76×` pure math); **`4.70 ms`** (`6.11×`) w/o `rceff`; **`N_FAIL = 0`** |
| **Kokkos-MONO** *(w/ `rceff`)* | **Pure Kokkos 4.7** | **+ $O(K)$ Monotone Proj.** | **`5.30 ± 0.04`** *(H200)* | **`5.42×`** (`+0.36 ms`); 100-step water **`52.7% → 9.28×10⁻¹⁴`** |
| **Native CUDA `v4.1` Mode 0→1** | CUDA C++17 | Original NSW6 (BB → ENG) | `31.91 → 11.70` *(H200)* | **`2.73×`** (Hopper eng); **`N_FAIL = 0`** |
| **Native CUDA `v4.1` Mode 2** | **CUDA C++17 + Graph** | **+ Our 3 Math Theorems** | **`4.07 ± 0.01`** *(H200)* | **`7.85×`** (`2.88×` math); **`3.89 ms`** (`8.21×`) w/o `rceff`; `85%` HBM3e |
| **Native CUDA `v4.1` Mode 3** | **CUDA C++17 + Graph** | **+ $O(K)$ Monotone Proj.** | **`4.41 ± 0.002`** *(H200)* | **`7.23×`** (`+0.35 ms`); 100-step water **`52.7% → 1.83×10⁻¹⁵`** |

> **Key takeaway**: Pure mathematical speedup is nearly identical in both frameworks (**`2.76×`** in Kokkos vs. **`2.88×`** in Native CUDA), and **100% pure portable Kokkos (`4.94 ms`) trails hand-tuned Hopper-only CUDA (`4.07 ms`) by only `0.88 ms` (`+21.6%`)**.

### Quick Start
```bash
cd baseline && shasum -a 256 -c SHA256SUMS                                    # Baseline integrity
cd ../work_opt/hopper_cuda && nvcc -O3 -std=c++17 -arch=sm_90a --prec-div=true --prec-sqrt=true bench_paper_ci.cu -o bench && ./bench
cd ../kokkos_bench && cmake -B build -DCMAKE_CXX_COMPILER=/path/to/nvcc_wrapper -DCMAKE_CXX_FLAGS="-O3 --prec-div=true --prec-sqrt=true" -DKokkos_ROOT=/path/to/Kokkos && cmake --build build -j && ./build/bench_kokkos_full
```

---

<a id="日本語"></a>
## 日本語

### 概要・研究目標・Todo
- **実装フレームワーク**：**ネイティブ CUDA C++17**（`sm_90a`・CUDA グラフ）および **100% 純ポータブル Kokkos 4.7**（`View` / `RangePolicy` / `TeamPolicy` 共有メモリ圧縮、生 CUDA 記述ゼロ）の両方を実装。数理改革は共通の Header-only C++ コア（`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`）に集約し、全カーネルで厳密倍精度 IEEE FP64・**レジスタスピル 0 B** を達成。
- **目標（Goal）**：**Ren らの論文（[`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313)）をガイドライン**として、NICAM NSW6 の **OpenACC / Triton / Kokkos / ネイティブ CUDA C++** を横断比較します。
- **Todo（`@QIAN JIAWEI`）**：
  - [x] **独占空闲 NVIDIA H200 での Kokkos 検証完了**：Ren 氏公開リポジトリ [`Masahiro-Ren/nicam_kokkos_microphysics`](https://github.com/Masahiro-Ren/nicam_kokkos_microphysics) を実測（1 領域 `16,641×96` 壁時計 `71.05 ms` / 純 GPU 核 `20.28 ms`）し、全領域（`166,410×80`）純 Kokkos 4.7 にて 8 ペアブロック評価を完了（**`28.71 → 13.64 → 4.94 ms`**、総合 **`5.81×`** / 純数理 **`2.76×`**、`N_FAIL = 0`）。
  - [ ] Miyabi-G（GH200）復旧後、**OpenACC / Triton / Kokkos / ネイティブ CUDA C++** の同一実チェックポイント比較を完了させる。

### 数理手法と効果
1. **雪モーメント（二進格子 $n \in \frac{1}{8}\mathbb{Z}$）＋乾燥格子選別**：桁落ち回避 `expm1` 倍増で 16 超越関数を `1 log + 2 exp` へ圧縮（`pow 131.5×`・`exp 7.31×`・`log 5.73×` 削減）＋乾燥格子 `83.3%` 除外（微物理処理格子 `6.0×` 削減）。
2. **鉛直沈降（フビニ因数分解＋$O(K)$ 単調包絡射影）**：共有質量重み $W_{k,j}=\sum_s L_s(k,j)a_{j,s}$ で 28 走査を集約（降水柱 `40%` のみ計算）。後向包絡射影 $\tilde{\xi}_k=\min(\xi_k,\tilde{\xi}_{k+1})$（`+0.35–0.36 ms`）で軌跡交差を全廃し、100 歩水質量誤差を **`52.7%` → 機械精度** へ低減。
3. **飽和調節（陪ラゲール 4 次ハウスホルダー法）**：厳密導関数 $q_{\mathrm{sat}}^{(m)}(T)=(-1)^m m!\,T^{-m} L_m^{(-\gamma)}(\alpha/T)\,q_{\mathrm{sat}}(T)$ でヤコビアン欠落（約 `5%`）を修復（最大 4 反復・反復数 `2.46×` 削減）。

### 独占空闲 H200 同一セッション実測比較（全領域 `166,410 × 80` = 1,331 万格子点、`dt = 30 s`、8 ペアブロック 95% CI）

| 実装 / モード | フレームワーク | アルゴリズム | 実行時間 (ms) | 速度比・検証 |
| :--- | :--- | :--- | :---: | :--- |
| **元ベースライン** | Fortran + OpenACC | 原 NSW6 | `78.19` *(GH200)* | 既存 GH200 ベースライン |
| **Triton 版** | Python / Triton | 原 NSW6 | `68.57` *(GH200)* | 既存 GH200 ベースライン |
| **Ren 氏 Kokkos（公開原本）** | Kokkos 4.7 | 原 NSW6（1 領域 `16,641×96`） | `71.05` 壁時計 / `20.28` GPU *(H200)* | 889 起動＋228 alloc/free（10 領域換算 GPU `~203 ms`） |
| **Kokkos-BB** | 純 Kokkos 4.7 | 原 NSW6（非融合・alloc なし） | `28.71 ± 0.09` *(H200)* | `1.00×` Kokkos 基準 |
| **Kokkos-ENG** | 純 Kokkos 4.7 | 原 NSW6（全領域融合） | `13.64 ± 0.14` *(H200)* | **`2.11×`**（実装最適化）・**`N_FAIL = 0`** |
| **Kokkos-MATH** *(`rceff` 含)* | **純 Kokkos 4.7** | **＋三大数理改革** | **`4.94 ± 0.02`** *(H200)* | **`5.81×`**（純数理 **`2.76×`**）、`rceff` 除 **`4.70 ms`**（`6.11×`）、**`N_FAIL = 0`** |
| **Kokkos-MONO** *(`rceff` 含)* | **純 Kokkos 4.7** | **＋$O(K)$ 単調射影** | **`5.30 ± 0.04`** *(H200)* | **`5.42×`**（`+0.36 ms`）、100 歩水誤差 **`52.7% → 9.28×10⁻¹⁴`** |
| **ネイティブ CUDA `v4.1` Mode 0→1** | CUDA C++17 | 原 NSW6（BB → ENG） | `31.91 → 11.70` *(H200)* | **`2.73×`**（Hopper 実装）・**`N_FAIL = 0`** |
| **ネイティブ CUDA `v4.1` Mode 2** | **CUDA C++17 + グラフ** | **＋三大数理改革** | **`4.07 ± 0.01`** *(H200)* | **`7.85×`**（純数理 **`2.88×`**）、`rceff` 除 **`3.89 ms`**（`8.21×`）、帯域 `85%` |
| **ネイティブ CUDA `v4.1` Mode 3** | **CUDA C++17 + グラフ** | **＋$O(K)$ 単調射影** | **`4.41 ± 0.002`** *(H200)* | **`7.23×`**（`+0.35 ms`）、100 歩水誤差 **`52.7% → 1.83×10⁻¹⁵`** |

> **要点**：純数理改革による高速化倍率は **Kokkos（`2.76×`）とネイティブ CUDA（`2.88×`）でほぼ一致**しており、**100% 純ポータブル Kokkos（`4.94 ms`）と Hopper 特化手書き CUDA（`4.07 ms`）の差はわずか `0.88 ms`（`+21.6%`）** に留まります。
