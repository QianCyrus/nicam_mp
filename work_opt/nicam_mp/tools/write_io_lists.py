#!/usr/bin/env python3
"""Emit INPUTS.txt / OUTPUTS.txt for each dataset package.

The filename alone does NOT tell you whether an array is an input or an
output: 19 files carry no _in/_out suffix and are inputs (metrics, aerosol,
scalars), while 14 others carry no suffix and are outputs (diagnostics with
no input counterpart, e.g. mpF_drv_precip, mpF_drv_cldfrac). One file even
looks like an input by prefix but is one by role only (mpF_pt_Vt). The
authority is the `role` field in manifest.json; these lists are that field
flattened into something a shell script or a Fortran driver can read
directly:

    grep -v '^#' INPUTS.txt      # -> one filename per line

Re-runnable:  python3 tools/write_io_lists.py
"""
import json
import os

ROOT = "/data/work/kokkos-nicam/nicam_mp/dataset"
INPUT_ROLES = ("metric", "input", "input-derived", "scalar")

# entry point -> (section title, filename prefix test)
SECTIONS = [
    ("static (read once; identical every step)",
     lambda n: not n.startswith(("mpF_drv_", "mpF_pt_", "mpF_sat_"))),
    ("mp_driver — the whole microphysics step (the mini-app's entry point)",
     lambda n: n.startswith("mpF_drv_")),
    ("precip_transport — sub-stage, only if you run it in isolation",
     lambda n: n.startswith("mpF_pt_")),
    ("saturation_adjustment — sub-stage, only if you run it in isolation",
     lambda n: n.startswith("mpF_sat_")),
]


def load_entries(manifest_path):
    m = json.load(open(manifest_path))
    arrs = m if isinstance(m, list) else m.get("arrays", m)
    if isinstance(arrs, dict):
        arrs = [dict(name=k, **v) for k, v in arrs.items()]
    return [(a.get("name"), a.get("role", "?")) for a in arrs]


def write_list(path, entries, want_input, ext, header):
    lines = [f"# {header}", "#",
             "# Roles come from manifest.json — the filename suffix is NOT a",
             "# reliable indicator. Parse with:  grep -v '^#' <this file>", "#"]
    total = 0
    for title, test in SECTIONS:
        names = sorted(n for n, r in entries
                       if test(n) and ((r in INPUT_ROLES) == want_input))
        if not names:
            continue
        lines.append(f"\n# --- {title}")
        lines += [n + ext for n in names]
        total += len(names)
    lines.insert(1, f"# {total} files")
    open(path, "w").write("\n".join(lines) + "\n")
    return total


for size in ("full", "region0"):
    for sub, ext in (("", ".bin"), ("npy", ".npy")):
        d = os.path.join(ROOT, size, sub)
        mf = os.path.join(d, "manifest.json")
        if not os.path.isfile(mf):
            print(f"skip (no manifest): {d}")
            continue
        e = load_entries(mf)
        ni = write_list(os.path.join(d, "INPUTS.txt"), e, True, ext,
                        f"INPUTS — cloud-microphysics mini-app dataset ({size})")
        no = write_list(os.path.join(d, "OUTPUTS.txt"), e, False, ext,
                        f"OUTPUTS (expected results for validation) ({size})")
        print(f"{size}/{sub or 'bin':7s}  INPUTS.txt {ni:3d}   OUTPUTS.txt {no:3d}")
