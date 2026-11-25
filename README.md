# NICAM NSW6 Cloud Microphysics — Mathematical Reformulation on Native CUDA & Portable Kokkos

<p align="center">
  <strong><a href="#english">English</a></strong> &nbsp;|&nbsp; <strong><a href="#日本語">日本語</a></strong>
</p>

<p align="center">
  <img alt="Frameworks: Native CUDA C++17 & Kokkos 4.7" src="https://img.shields.io/badge/Frameworks-Native%20CUDA%20C%2B%2B17%20%7C%20Kokkos%204.7-76B900?style=flat-square" />
  <img alt="Precision: Strict IEEE FP64" src="https://img.shields.io/badge/Precision-Strict%20IEEE%20FP64-00599C?style=flat-square" />
  <img alt="Spill: 0 Bytes" src="https://img.shields.io/badge/Register%20Spill-0%20Bytes-2EA44F?style=flat-square" />
  <img alt="Oracle: N_FAIL = 0" src="https://img.shields.io/badge/Validation-N__FAIL%20%3D%200-2EA44F?style=flat-square" />
  <img alt="Speedup: 7.86x CUDA / 5.14x Kokkos" src="https://img.shields.io/badge/H200%20Speedup-7.86%C3%97%20CUDA%20%7C%205.14%C3%97%20Kokkos-E36209?style=flat-square" />
</p>

---

<a id="english"></a>
## English

### Overview, Goal & Todo
- **Frameworks**: Implemented in both **Native NVIDIA CUDA C++17** (`sm_90a`, CUDA Graphs) and **100% Pure Portable Kokkos 4.7** (`View` / `RangePolicy` / `TeamPolicy` scratch compaction, zero raw CUDA syntax), sharing one header-only C++ math core (`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`), all in strict IEEE FP64 (`--prec-div=true --prec-sqrt=true`, no `--use_fast_math`).
- **Goal**: Following **Ren et al. ([`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313))** as a guideline, compare **Fortran/OpenACC**, **Triton**, **C++/Kokkos**, and **Native CUDA C++** for NICAM NSW6.
- **Todo (`@QIAN JIAWEI`)**:
  - [x] Audit & benchmark **Kokkos** on NVIDIA H200: profiled upstream [`Masahiro-Ren/nicam_kokkos_microphysics`](https://github.com/Masahiro-Ren/nicam_kokkos_microphysics) (`21.60 ms` pure GPU kernels for 1 region `16,641×96`, 889 launches) and validated our mathematical reductions inside **100% pure portable Kokkos 4.7** (`N_FAIL = 0`, 0 B spill, **`5.14×`** pure-math GPU kernel speedup).
  - [ ] Run unified **Miyabi-G (GH200)** real-checkpoint comparisons across **OpenACC / Triton / Kokkos / Native CUDA C++** post-maintenance.

### Mathematical Methods & Per-Stage Speedup (Same Kokkos 4.7 Framework, H200 `166,410×80`)
1. **Snow Moments (Dyadic Lattice $n \in \frac{1}{8}\mathbb{Z}$) + Dry Screening**: Cancellation-free `expm1`-doubling reduces 16 transcendentals (`8 exp + 8 pow`) per cell to `1 log + 2 exp` (`131.5×` fewer `pow`, `7.31×` fewer `exp`), plus `83.3%` dry-cell screening → **Microphysics core: `15.12 → 1.71 ms` (`8.86×`)**.
2. **Sedimentation (Fubini Factorisation + $O(K)$ Monotone Projection)**: Shared mass weight $W_{k,j}=\sum_s L_s(k,j)a_{j,s}$ collapses 28 flux sweeps + `40%` precipitating-column compaction → **Full-level: `5.65 → 1.25 ms` (`4.52×`), Half-level: `1.31 → 0.40 ms` (`3.30×`)**; backward-envelope projection $\tilde{\xi}_k=\min(\xi_k,\tilde{\xi}_{k+1})$ removes all trajectory crossings (`+0.34 ms`).
3. **Saturation Adjustment (Associated-Laguerre Householder-4)**: Exact derivatives $q_{\mathrm{sat}}^{(m)}(T)=(-1)^m m!\,T^{-m} L_m^{(-\gamma)}(\alpha/T)\,q_{\mathrm{sat}}(T)$ fix a `~5%` Jacobian defect (`≤ 4` iters, `2.46×` fewer) → **SatAdj+Diag: `3.56 → 1.09 ms` (`3.26×`)**.

### Cross-Framework Performance & Conservation (`166,410 × 80` = 13.31M cells, `dt = 30 s`)

| Implementation / Mode | Framework | Algorithm | Step Time | Notes / Validation |
| :--- | :--- | :--- | :---: | :--- |
| **Baseline** | Fortran + OpenACC | Original NSW6 | `78.19 ms` *(GH200)* | Recorded baseline |
| **Triton** | Python / Triton | Original NSW6 | `68.57 ms` *(GH200)* | Recorded baseline |
| **Ren Kokkos (upstream)** | Kokkos 4.7 | Original NSW6 (`16,641×96`, 1 reg) | `21.60 ms` *(H200 GPU sum)* | 889 launches (`~180–216 ms` 10-reg eq.) |
| **Kokkos-BB / ENG** | Pure Kokkos 4.7 | Original NSW6 (whole-domain) | `26.19 ms` *(unpreempted)* / `31.32 ms` *(shared)* | **`N_FAIL = 0`**; `Kokkos-BB`: `65.41 ms` shared |
| **Kokkos-MATH** *(w/ `rceff`)* | **Pure Kokkos 4.7** | **+ Our 3 Math Theorems** | **`~5.10 ms`** *(unpreempted)* / `12.56 ms` *(shared)* | **`5.14×`** math GPU (`5.21×` vs BB); **`N_FAIL = 0`** |
| **Kokkos-MONO** *(w/ `rceff`)* | **Pure Kokkos 4.7** | **+ $O(K)$ Monotone Proj.** | **`~5.44 ms`** *(+0.34 ms)* / `12.90 ms` *(shared)* | 100-step water drift **`52.7% → 9.28×10⁻¹⁴`** |
| **Native CUDA `v4.1` Mode 2** | **CUDA C++17 + Graph** | **+ Our 3 Math Theorems** | **`4.06 ms`** *(H200 excl.)* / `3.89 ms` *(no `rceff`)* | **`7.86×`** vs BB-H200 (`31.90 ms`); `85%` HBM3e |
| **Native CUDA `v4.1` Mode 3** | **CUDA C++17 + Graph** | **+ $O(K)$ Monotone Proj.** | **`4.41 ms`** *(H200 excl., +0.36 ms)* | **`7.23×`**; 100-step water **`52.7% → 1.83×10⁻¹⁵`** |

> *Measurement note*: Kokkos H200 tests ran on a shared H200 under concurrent tenant load (`690 W / 100% util`), so both paired shared-load wall-clock (`8×30` blocks) and unpreempted `nsys Min` GPU kernel sums are reported transparently. Pure portable Kokkos (`~5.10 ms`) trails hand-tuned Hopper-only CUDA (`4.06 ms`) by only `~25%`.

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
- **実装フレームワーク**：**ネイティブ CUDA C++17**（`sm_90a`・CUDA グラフ）および **100% 純ポータブル Kokkos 4.7**（`View` / `RangePolicy` / `TeamPolicy` 共有メモリ圧縮、生 CUDA 記述ゼロ）の両方で実装。3 つの数理改革は共通の Header-only C++ コア（`work_opt/nicam_mp/cuda_kernels/nsw6_math_core.cuh`）に集約し、すべて厳密倍精度 IEEE FP64（`--prec-div=true --prec-sqrt=true`）で検証済みです。
- **目標（Goal）**：**Ren らの論文（[`arXiv:2610.41313`](https://arxiv.org/abs/2610.41313)）をガイドライン**として、NICAM NSW6 の **OpenACC / Triton / Kokkos / ネイティブ CUDA C++** を横断比較します。
- **Todo（`@QIAN JIAWEI`）**：
  - [x] **H200 での Kokkos 版性能検証完了**：Ren 氏公開リポジトリ [`Masahiro-Ren/nicam_kokkos_microphysics`](https://github.com/Masahiro-Ren/nicam_kokkos_microphysics) を解析（1 領域 `16,641×96` 純 GPU 核合計 `21.60 ms`・889 起動）し、**純 Kokkos 4.7 同フレームワーク内**に我々の数理手法を適用して **`N_FAIL = 0`**・スピル 0 B・**純数理 `5.14 倍`** 高速化を実証。
  - [ ] Miyabi-G（GH200）復旧後、**OpenACC / Triton / Kokkos / ネイティブ CUDA C++** の同一実チェックポイント比較を完了させる。

### 数理手法とステージ別高速化（同一 Kokkos 4.7 フレームワーク、H200 `166,410×80`）
1. **雪モーメント（二進格子 $n \in \frac{1}{8}\mathbb{Z}$）＋乾燥格子選別**：桁落ち回避 `expm1` 倍増で 1 格子点 16 回の超越関数を `1 log + 2 exp` に圧縮（`pow 131.5×`・`exp 7.31×` 削減）＋乾燥格子 `83.3%` 除外 → **微物理コア：`15.12 → 1.71 ms`（`8.86×`）**。
2. **鉛直沈降（フビニ因数分解＋$O(K)$ 単調包絡射影）**：共有質量重み $W_{k,j}=\sum_s L_s(k,j)a_{j,s}$ で 28 系統走査を集約＋降水柱 `40%` 圧縮 → **全層沈降：`5.65 → 1.25 ms`（`4.52×`）、半層沈降：`1.31 → 0.40 ms`（`3.30×`）**。後向包絡射影 $\tilde{\xi}_k=\min(\xi_k,\tilde{\xi}_{k+1})$ で軌跡交差を完全除去（`+0.34 ms`）。
3. **飽和調節（陪ラゲール 4 次ハウスホルダー法）**：厳密導関数 $q_{\mathrm{sat}}^{(m)}(T)=(-1)^m m!\,T^{-m} L_m^{(-\gamma)}(\alpha/T)\,q_{\mathrm{sat}}(T)$ でヤコビアン欠落（約 `5%`）を修復（最大 4 反復・`2.46×` 削減） → **飽和調節＋診断：`3.56 → 1.09 ms`（`3.26×`）**。

### フレームワーク横断性能・保存性サマリー（全領域 `166,410 × 80` = 1,331 万格子点、`dt = 30 s`）

| 実装 / モード | フレームワーク | アルゴリズム | 1 ステップ実行時間 | 備考・検証 |
| :--- | :--- | :--- | :---: | :--- |
| **元ベースライン** | Fortran + OpenACC | 原 NSW6 | `78.19 ms` *(GH200)* | 計測済ベースライン |
| **Triton 版** | Python / Triton | 原 NSW6 | `68.57 ms` *(GH200)* | 計測済ベースライン |
| **Ren 氏 Kokkos（公開原本）** | Kokkos 4.7 | 原 NSW6（1 領域 `16,641×96`） | `21.60 ms` *(H200 核合計)* | 889 起動（10 領域換算 `~180–216 ms`） |
| **Kokkos-BB / ENG** | 純 Kokkos 4.7 | 原 NSW6（全領域） | `26.19 ms` *(無競合核)* / `31.32 ms` *(共有負荷)* | **`N_FAIL = 0`**、`Kokkos-BB`: `65.41 ms` |
| **Kokkos-MATH** *(`rceff` 含)* | **純 Kokkos 4.7** | **＋三大数理改革** | **`~5.10 ms`** *(無競合核)* / `12.56 ms` *(共有負荷)* | **純数理 `5.14×`**（対 BB `5.21×`）・**`N_FAIL = 0`** |
| **Kokkos-MONO** *(`rceff` 含)* | **純 Kokkos 4.7** | **＋$O(K)$ 単調射影** | **`~5.44 ms`** *(+0.34 ms)* / `12.90 ms` *(共有負荷)* | 100 歩水質量誤差 **`52.7% → 9.28×10⁻¹⁴`** |
| **ネイティブ CUDA `v4.1` Mode 2** | **CUDA C++17 + グラフ** | **＋三大数理改革** | **`4.06 ms`** *(H200 占有)* / `3.89 ms` *(`rceff` 除)* | 対 BB-H200 **`7.86×`**（実効帯域 `85%`） |
| **ネイティブ CUDA `v4.1` Mode 3** | **CUDA C++17 + グラフ** | **＋$O(K)$ 単調射影** | **`4.41 ms`** *(H200 占有, +0.36 ms)* | **`7.23×`**・100 歩水誤差 **`52.7% → 1.83×10⁻¹⁵`** |

> *注記*：H200 の Kokkos 計測時は他ジョブとの共有稼働下（`690 W / 100%`）であったため、同負荷ペア計測値と `nsys Min` 無競合核合計の両方を併記しています。純ポータブル Kokkos（`~5.10 ms`）は Hopper 特化の手書き CUDA（`4.06 ms`）と約 `25%` 差に迫り、数理改革の効果がフレームワーク非依存であることを実証しました。
