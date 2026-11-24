# NICAM Cloud-Microphysics Mini-App

A standalone, single-step benchmark and validation harness for NICAM's NSW6
cloud microphysics, extracted from the full model.

It loads a real atmospheric state (dumped from a production NICAM run), calls
the microphysics driver `mp_driver` exactly once — condensation, freezing,
melting, collection, sedimentation, saturation adjustment, 30 model seconds of
cloud physics — and compares every output array against the recorded reference.
No MPI decomposition, no NetCDF, no icosahedral grid machinery: the
microphysics is column-local, so to this program the atmosphere is just
"N columns × 80 levels".

Use it to:
- **port** the NSW6 scheme to another language/framework (JAX, CUDA, ...) with
  a bit-level reference to validate against, or
- **benchmark** the Fortran OpenACC implementation on a GPU (166k columns of
  genuinely convecting atmosphere — 44 % of columns are precipitating).

The sources are copied from NICAM with original headers and comments kept;
every deviation is marked with a greppable `!--- MINIAPP:` comment
(`grep -rn MINIAPP nhm/`).

# Package Layout

```
nicam_mp/
├── nhm/
│   ├── prg_mp_miniapp.f90     driver: load → mp_driver → validate → report
│   ├── physics/               the ported kernel (5 files, NICAM originals)
│   │   ├── mod_mp_driver.f90        entry point / process orchestration
│   │   ├── mod_mp_nsw6.f90          NSW6 6-category cold-rain core
│   │   ├── mod_precip_transport.f90 sedimentation / precipitation flux
│   │   ├── mod_vadv1d.f90           semi-Lagrangian vertical advection
│   │   └── mod_satadjust.f90        saturation adjustment (Newton)
│   └── share/                 support modules (constants, thermodynamics, ...)
├── dataset/
│   ├── region0/               16,641 columns, 2.2 GB  — develop & debug
│   ├── full/                  166,410 columns, 12 GB  — performance runs
│   └── vgrid78.dat            vertical grid definition
├── run/
│   ├── local/                 docker runner + configs (this machine)
│   └── miyabi-g/              PBS job script + configs (Miyabi-G)
├── docker/Dockerfile          reproduces the local build image
├── compile_docker.sh          local GPU build   (inside docker)
├── compile_docker_cpu.sh      local CPU build   (inside docker)
├── compile.sh                 Miyabi-G GPU build (login node, make -j 8)
├── compile_cpu.sh             Miyabi-G CPU build
├── Makefile / Mkinclude / sysdep/   NICAM-style build system
├── tools/                     dataset packaging / npy conversion scripts
└── README.md                  this file
```

Binaries are installed to `./bin` and copied to `./bin_cpu` / `./bin_gpu` by
the compile scripts.

# Dataset

| package   | columns  | size   | intended use |
|-----------|----------|--------|--------------|
| `dataset/region0` | 16,641  | 2.2 GB | development, debugging, CI — GPU ~1 s, CPU ~10 s |
| `dataset/full`    | 166,410 | 12 GB  | performance evaluation — saturates a modern GPU |

Two file formats, same values:

- `*.bin` — big-endian Fortran unformatted records. **This is what the Fortran
  mini-app reads.**
- `npy/*.npy` — numpy mirror (float64, little-endian, Fortran order). **This is
  the input for Python/JAX ports** — no Fortran record parsing needed.

Each package also carries `manifest.json` (per array: shape, units, meaning,
statistics) and `INPUTS.txt` / `OUTPUTS.txt` — the authoritative lists of what
is an input and what is an expected output. **Do not guess roles from the
filenames** (some outputs have no `_in`/`_out` suffix); see
"Input and Output List" at the end of this document.

**Where the data comes from**: a production NICAM run (gl07rl00, full physics)
was integrated for 10 model days; the restart was then advanced by exactly one
dynamics step and handed to `mp_driver` — i.e. the recorded input/output pair
is the model's own microphysics call at production step 28801, on an
atmosphere with active convection. `region0` is the first of the ten regions
of `full`. The reference was produced by the **CPU** build, which is why a CPU
port can (and must) match it bit-for-bit. Per-array documentation (shape, units,
meaning, statistics): `manifest.json` inside each dataset package.

# How to build

## docker and dockerfile

Local builds run inside the `fortran-nvhpc:26.3` container (NVHPC 26.3
compilers + NetCDF stack). If you don't have the image:

```bash
docker build -t fortran-nvhpc:26.3 docker/
```

(Builds HDF5 + NetCDF-C/Fortran from source with nvc/nvfortran — takes a
while, needed once.)

## how to build on local

```bash
./compile_docker_cpu.sh     # CPU/OpenMP build  -> ./bin_cpu/mp_miniapp
./compile_docker.sh         # GPU/OpenACC build -> ./bin_gpu/mp_miniapp
```

A `make allclean` runs automatically between CPU↔GPU builds (object flags
differ; the scripts handle it).

## how to build on miyabi-g

On a login node (the scripts cap `make -j 8` — unbounded parallel builds get
killed by the login-node memory limit, with no error message):

```bash
./compile.sh          # GPU (OpenACC, cc90)  -> ./bin_gpu/mp_miniapp
./compile_cpu.sh      # CPU (OpenMP)         -> ./bin_cpu/mp_miniapp
```

Both load `nvidia/25.9` + `nv-hpcx` and use `sysdep/Makedef.Miyabi-G`.

# How to run

## how to run on local

```bash
cd run/local
./run_docker.sh                                            # CPU, region0
GPU=T BIN=bin_gpu CNF=mp_miniapp_gpu.cnf  ./run_docker.sh  # GPU, region0
GPU=T BIN=bin_gpu CNF=mp_miniapp_full.cnf ./run_docker.sh  # GPU, full size
```

Env knobs: `BIN` (binary dir), `CNF` (config), `GPU=T` (pass the GPU to
docker), `OMP` (threads). The log is written to `msg_mp.pe000000`
(+ binary results in `output/` if `write_output` is enabled).

## how to run on miyabi-g

```bash
cd run/miyabi-g
# edit run-MYBG.sh once: set your  #PBS -W group_list
qsub run-MYBG.sh                   # GPU, region0, GPU-tolerance config
```

Variants: `qsub -v BIN=../../bin_cpu,CNF=mp_miniapp.cnf run-MYBG.sh` for the
bit-exact CPU check; `CNF=mp_miniapp_full.cnf` for the full-size performance
run (raise the walltime).

# How to read the msg.pe000000

## The validation table

```
                   ARRAY STATUS    MAX|DIFF|      MAX_REL      N_FAIL        N_ELEM
rhog                     OK       0.0000E+00   0.0000E+00           0       1331280
q:qc                     OK       9.9818E-12   2.2253E-02           0       1331280
```

| column | meaning |
|---|---|
| `ARRAY`  | variable name as NICAM's `mod_physics.f90` names it at the `mp_driver` call site. Multi-category arrays are split per water species: `q:qv, q:qc, q:qr, q:qi, q:qs, q:qg` (vapour, cloud, rain, ice, snow, graupel) and `precip_mp:liq / :sol` |
| `STATUS` | `OK` / `DIFFER` — verdict at the configured tolerance |
| `MAX\|DIFF\|` | largest absolute difference vs the reference |
| `MAX_REL` | largest relative difference (references below the model epsilon 1e-16 are excluded from the division) |
| `N_FAIL` | number of elements breaching the tolerance — the triage number: a handful ⇒ isolated edge cells; ~half the field ⇒ real bug |
| `N_ELEM` | number of elements compared |

**What to expect.** A CPU float64 port compiled like the reference scores
**exactly 0.0 everywhere** (config `mp_miniapp.cnf`: rel 1e-13 / abs 0). A GPU
build sits at the compiler fast-math floor (config `mp_miniapp_gpu.cnf`:
rel 1e-9 / abs 1e-10): typically `pre` ~2e-11 relative, and `q:qc` may show
`MAX_REL` up to ~2e-2 on cloud-edge cells whose reference value is a tiny
cancellation residual — while its absolute difference stays ~1e-11 kg/kg.
That is compiler arithmetic, not an error (the reference compiler moves more
than that against itself when its strict-IEEE flag is toggled). If a run
FAILs the GPU gate, something real is wrong — for JAX ports the classic cause
is float32: set `jax.config.update("jax_enable_x64", True)` before creating
any array.

## The time report

NICAM-style rap timers at the end of the log:

```
*** ID=002 : INIT_Load_input        time to read the dataset (I/O — exclude from perf numbers)
*** ID=003 : MAIN_Main_Loop         everything below
*** ID=004 : MAIN_Restore_state     state reset between repeats (bookkeeping)
*** ID=005 : MAIN_MP_driver         THE number: microphysics compute time
*** ID=006 : MAIN_____Precip_Transport        sub-timer inside mp_driver
*** ID=007 : MAIN_____Saturation_Adjustment   sub-timer inside mp_driver
*** ID=008 : CHCK_Validate          comparison against the reference
```

Per-step cost = `MAIN_MP_driver` T(avg) ÷ `num_repeat`. Reference numbers on
the full dataset, one device: ~10.9 s/step CPU (1 thread) vs ~0.16 s/step GPU
(RTX PRO 6000 Blackwell), ≈ 68×.

**GPU timing**: the `MAIN_MP_driver` total is always exact. For real GPU time
in the sub-section rows (`Precip_Transport`, `Saturation_Adjustment`), set
`PROF_rap_sync = .true.` in `&PROFPARAM`.

# The configuration file

`run/*/mp_miniapp*.cnf`, NICAM namelist format. What a user may modify:

- `&MP_MINIAPPPARAM` — the mini-app itself:
  - `dataset_dir`, `output_dir` — dataset location (relative to the run dir)
  - `num_repeat` — how many times to re-run the same step (use ≥3 for timing;
    the state is restored between repeats, so every iteration does identical
    work)
  - `tolerance`, `tolerance_abs` — the validation gate. **Use the shipped
    pairs, don't invent numbers**: `mp_miniapp.cnf` = bit-exact CPU gate
    (1e-13 / 0), `mp_miniapp_gpu.cnf` = measured GPU envelope (1e-9 / 1e-10,
    rationale in the file comments)
  - `validate`, `write_output` — toggle the comparison / binary output
- `&ADMPARAM lall` — 1 for `region0`, 10 for `full`; must match `dataset_dir`
- `&PROFPARAM PROF_rap_level` — timer detail (2 = the table above);
  `PROF_rap_sync = .true.` — real GPU time in the per-section rows
- Everything else (`&RUNCONFPARAM`, `&NM_MP_NSW6`, `&GRDPARAM`) reproduces the
  production physics configuration — the reference outputs were computed with
  exactly these values, so change them only if you know why.

# Simple call graph

```
prg_mp_miniapp                          (driver — plays the role of NICAM's physics_step)
 ├── load inputs                        31 arrays (see below)
 └── for each repeat:
      restore pristine state
      for l = 1, lall:                  (column blocks, independent)
        mp_driver                       mod_mp_driver.f90 — the entry point
         ├── negative_filter            remove negative water left by advection
         ├── mp_nsw6                    mod_mp_nsw6.f90 — all microphysical
         │                              conversion rates (49 processes) + Vt
         ├── precip_transport_new       mod_precip_transport.f90 — sedimentation
         │    └── vadv1d_getflux_new    mod_vadv1d.f90 — semi-Lagrangian fall
         └── SATURATION_adjustment      mod_satadjust.f90 — Newton iteration,
                                        restores liquid/vapour equilibrium
      validate                          compare the 23 outputs
      report                            table + rap timers
```

Porting note: the two data-dependent constructs that don't map 1:1 onto array
frameworks are the Newton loop in `SATURATION_adjustment` (iterate-until-
converged with block-exit semantics) and the departure-cell search in
`vadv1d` — both need masked fixed-iteration forms (e.g. `lax.while_loop`).

# Input and Output List

Machine-readable versions live in every package directory: `INPUTS.txt` and
`OUTPUTS.txt` (grep out the `#` lines to get one filename per row); shapes,
units and statistics per array are in `manifest.json`.

## Inputs (31 for a whole-step run)

**Static — read once, identical every step (16):**

| file | meaning |
|---|---|
| `mpF_z`, `mpF_zh` | height of full / half model levels |
| `mpF_phi` | geopotential |
| `mpF_gam2`, `mpF_gam2h`, `mpF_gsgam2`, `mpF_gsgam2h` | vertical metric factors |
| `mpF_ix mpF_iy mpF_iz mpF_jx mpF_jy mpF_jz` | local horizontal unit vectors |
| `mpF_unccn` | cloud-condensation nuclei (production default fill) |
| `mpF_aop_dnovrn` | aerosol optical parameter (production default fill) |
| `mpF_params` | scalars: `dt = 30 s`, current model time |

**Per-step state — what the microphysics changes (15, `mpF_drv_*_in`):**

| file(s) | meaning |
|---|---|
| `rhog` | density (×G½γ², the metric-weighted prognostics) |
| `rhogvx rhogvy rhogvz rhogw` | momentum |
| `rhoge` | internal energy |
| `rhogq` | tracer mass (6 water species) |
| `rho pre tem` | density, pressure, temperature (diagnostics) |
| `vx vy vz w` | velocities |
| `q` | tracer mixing ratios (6 water species) |

(`mpF_pt_*` / `mpF_sat_*` files in `region0` belong to the two sub-stage entry
points — only needed when debugging precip-transport or saturation-adjustment
in isolation; `full` does not carry them.)

## Outputs (23 expected results, validated)

| file(s) | meaning |
|---|---|
| `mpF_drv_*_out` (15) | the updated state — same variables as the 15 inputs |
| `mpF_drv_precip` | surface precipitation flux, liquid / solid [kg/m²/s] |
| `mpF_drv_precip_lh` | latent-heat flux carried by precipitation |
| `mpF_drv_precip_rhoein` | internal energy carried by precipitation |
| `mpF_drv_precip_rhophi` | potential energy carried by precipitation |
| `mpF_drv_precip_rhokin` | kinetic energy carried by precipitation |
| `mpF_drv_precip_energy` | total energy flux of precipitation |
| `mpF_drv_cldfrac` | diagnosed cloud fraction |
| `mpF_drv_qhyd` | total hydrometeor mixing ratio |
