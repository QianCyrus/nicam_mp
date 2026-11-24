#!/usr/bin/env python3
"""dump_to_npy.py -- package the kksref Stage-F microphysics dump.

Reads the raw Stage-F dump written by the Fortran harness
(nicam/NICAM/nhm/driver/prg_kksref.f90, do_stageF=.true.) and produces the
self-contained mini-app dataset:

    <out>/full/            all 10 regions   (production size)
    <out>/full/*.bin         raw Fortran sequential-unformatted, BIG-ENDIAN
    <out>/full/manifest.txt  inventory + roles + units + per-array stats
    <out>/full/npy/*.npy     numpy mirror (native little-endian float64)
    <out>/full/npy/manifest.json

    <out>/region0/         region l=1 only  (debug/CI size)
      ... same layout ...

Every Stage-F array has the region index l as its LAST (slowest) Fortran
dimension, so region 0 is simply the first 1/ADM_lall of the record payload --
no re-run of the harness is needed.

The raw .bin files are byte-faithful to what NICAM wrote (4-byte big-endian
record markers, column-major payload, IEEE double); for region0 the payload is
sliced and the markers are rewritten to the new length.

Usage:
    ./dump_to_npy.py [--dump run/mp-dataset/dump] [--out nicam_mp/dataset]
                     [--sizes full,region0] [--lall 10] [--no-raw] [--no-npy]
"""

import argparse
import json
import os
import struct
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mp_arrays import ARRAYS, SCALARS, PROVENANCE, DIMS  # noqa: E402

BE = {"f8": ">f8", "i4": ">i4"}
NATIVE = {"f8": "<f8", "i4": "<i4"}
PREFIX = "mpF_"


# --------------------------------------------------------------------------
def read_manifest(dirname):
    """Parse the harness manifest.txt -> {name: (dtype, dims)} (insertion order)."""
    path = os.path.join(dirname, "manifest.txt")
    if not os.path.isfile(path):
        sys.exit(f"ERROR: no manifest.txt in {dirname}")
    entries = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                continue
            name, dtype, rank = parts[0], parts[1], int(parts[2])
            entries[name] = (dtype, tuple(int(x) for x in parts[3:3 + rank]))
    return entries


def read_record(path, dtype, dims):
    """Read one Fortran sequential-unformatted big-endian record -> ndarray (F order)."""
    nbytes = int(np.prod(dims)) * np.dtype(BE[dtype]).itemsize
    with open(path, "rb") as fh:
        (marker,) = struct.unpack(">i", fh.read(4))
        if marker != nbytes:
            raise IOError(f"{path}: record marker {marker} != expected {nbytes}")
        data = np.frombuffer(fh.read(marker), dtype=BE[dtype])
        (trailer,) = struct.unpack(">i", fh.read(4))
        if trailer != marker:
            raise IOError(f"{path}: trailing record marker mismatch")
    return data.reshape(dims, order="F")


def write_record(path, arr, dtype):
    """Write one Fortran sequential-unformatted big-endian record."""
    payload = np.asfortranarray(arr).astype(BE[dtype], copy=False)
    raw = payload.tobytes(order="F")
    with open(path, "wb") as fh:
        fh.write(struct.pack(">i", len(raw)))
        fh.write(raw)
        fh.write(struct.pack(">i", len(raw)))


def stats(arr):
    a = np.asarray(arr, dtype=np.float64).ravel()
    finite = np.isfinite(a)
    n = a.size
    return {
        "count": int(n),
        "min": float(a.min()) if n else 0.0,
        "max": float(a.max()) if n else 0.0,
        "mean": float(a.mean()) if n else 0.0,
        "nonzero_fraction": float(np.count_nonzero(a) / n) if n else 0.0,
        "nonfinite_count": int(n - int(finite.sum())),
    }


TRACERS = ["qv", "qc", "qr", "qi", "qs", "qg"]
PHASES = ["liquid", "solid"]


def slice_stats(arr, dims, dim_labels):
    """Per-tracer (v) or per-phase stats for arrays that carry such an axis."""
    if "v" in dim_labels:
        axis, labels = dim_labels.index("v"), TRACERS
    elif "phase" in dim_labels:
        axis, labels = dim_labels.index("phase"), PHASES
    else:
        return None
    out = {}
    for i, lab in enumerate(labels):
        out[lab] = stats(np.take(arr, i, axis=axis))
    return out


def meta_for(name):
    m = ARRAYS.get(name)
    if m is None:
        return {"role": "unknown", "units": "?", "desc": "(no metadata entry)"}
    return m


# --------------------------------------------------------------------------
def build(dump_dir, out_root, size, lall, do_raw, do_npy, exclude=()):
    entries = {k: v for k, v in read_manifest(dump_dir).items() if k.startswith(PREFIX)}
    if exclude:
        entries = {k: v for k, v in entries.items()
                   if not any(k.startswith(p) for p in exclude)}
    if not entries:
        sys.exit(f"ERROR: no {PREFIX}* arrays in {dump_dir}/manifest.txt "
                 "(was the harness run with do_stageF=.true.?)")

    size_dir = os.path.join(out_root, size)
    npy_dir = os.path.join(size_dir, "npy")
    os.makedirs(size_dir, exist_ok=True)
    if do_npy:
        os.makedirs(npy_dir, exist_ok=True)

    records = []
    for name, (dtype, dims) in entries.items():
        arr = read_record(os.path.join(dump_dir, name + ".bin"), dtype, dims)

        sliced = False
        if size == "region0" and dims and dims[-1] == lall:
            arr = np.asfortranarray(arr[..., 0:1])   # keep the l axis, length 1
            dims = tuple(dims[:-1]) + (1,)
            sliced = True

        if do_raw:
            write_record(os.path.join(size_dir, name + ".bin"), arr, dtype)
        if do_npy:
            np.save(os.path.join(npy_dir, name + ".npy"),
                    np.asfortranarray(arr.astype(NATIVE[dtype])))

        m = meta_for(name)
        labels = dim_names(name, dims, lall if not sliced else 1)
        records.append({
            "name": name,
            "role": m["role"],
            "shape": list(dims),
            "dim_names": labels,
            "dtype": "float64" if dtype == "f8" else "int32",
            "bin_byte_order": "big-endian",
            "npy_byte_order": "little-endian (numpy native)",
            "order": "Fortran (column-major); first index varies fastest",
            "units": m["units"],
            "description": m["desc"],
            "bin_bytes": int(np.prod(dims)) * 8 + 8,
            "stats": stats(arr),
            "slice_stats": slice_stats(arr, dims, labels),
        })
        del arr

    return size_dir, npy_dir, records


def dim_names(name, dims, lall):
    """Best-effort dimension labels for a Stage-F array shape."""
    if name == "mpF_params":
        return ["param"]
    names = []
    for i, d in enumerate(dims):
        last = (i == len(dims) - 1)
        if d == 16641:
            names.append("gall")
        elif d == 80:
            names.append("k")
        elif d == 6:
            names.append("v")
        elif last and (d == lall or d == 1):
            names.append("l")
        elif d == 2:
            names.append("phase")
        elif d == 1:
            names.append("sfc")
        else:
            names.append(f"n{d}")
    return names


# --------------------------------------------------------------------------
def write_manifests(size_dir, npy_dir, records, size, dump_dir, do_npy):
    generated = time.strftime("%Y-%m-%d %H:%M:%S %z")
    total_bin = sum(r["bin_bytes"] for r in records)

    doc = {
        "dataset": f"NICAM NSW6 cloud-microphysics mini-app dataset ({size})",
        "size_variant": size,
        "generated": generated,
        "source_dump": dump_dir,
        "provenance": PROVENANCE,
        "dimensions": DIMS,
        "formats": {
            "raw": "<name>.bin -- ONE Fortran sequential-unformatted record: "
                   "4-byte big-endian record marker (payload length in bytes), "
                   "payload, same 4-byte marker again. Payload is IEEE "
                   "binary64 BIG-ENDIAN in Fortran column-major order. This is "
                   "exactly what NICAM wrote (built with -Mbyteswapio); a "
                   "Fortran mini-app reads it with a plain "
                   "open(...,form='unformatted') + read(fid) arr.",
            "npy": "npy/<name>.npy -- numpy v1 .npy, float64, little-endian "
                   "(native), fortran_order=True, identical values and shape.",
        },
        "region_layout": "The region index l is the LAST Fortran dimension of "
                         "every array, so region N is the contiguous slice "
                         "[N*size/l_total : (N+1)*size/l_total] of the payload. "
                         "region0 was produced by taking that slice of the full "
                         "dump; no second harness run was needed.",
        "entry_points": {
            "1_mp_driver": {
                "call": "mp_driver(ijdim,kdim,l_region, rhog,rhogvx,rhogvy,rhogvz,"
                        "rhogw,rhoge,rhogq, rho,pre,tem,vx,vy,vz,w,q, precip,...,"
                        "GDCFRC,GDCLW, CBMFX,unccn,factor_n_aie,rhogqke, z,zh,"
                        "ix..jz, phi,gam2,gam2h,gsgam2,gsgam2h, dt,ct)",
                "in": "mpF_drv_*_in + all metric arrays + mpF_unccn + "
                      "mpF_aop_dnovrn (factor_n_aie) + mpF_params",
                "out": "mpF_drv_*_out + mpF_drv_precip* + mpF_drv_cldfrac + "
                       "mpF_drv_qhyd",
                "note": "One call per region = the whole MP step: negative-value "
                        "filter -> NSW6 core -> precipitation transport -> "
                        "saturation adjustment -> cloud-fraction diagnosis. "
                        "CBMFX and rhogqke enter as zero and are not dumped. "
                        "CBMFX=0 is production-exact (cp_driver('NONE') "
                        "zero-fills it, mod_cp_driver.f90:558). rhogqke is "
                        "nonzero in production (restart diagvar qke) but is a "
                        "DEAD input for MP_TYPE='NSW6': mp_driver reads it only "
                        "in the NDW6 branch (mod_mp_driver.f90:1217), and CBMFX "
                        "only in the LSC branch (:699) - neither affects any "
                        "NSW6 output.",
            },
            "2_precip_transport_new": {
                "call": "precip_transport_new(ijdim, rhog,...,q,qd, WPREC, "
                        "precipitating_flag, z,zh, precip,precip_rhoe,"
                        "precip_lh_heat,precip_rhophi,precip_rhokin,precip_flux, "
                        "gsgam2,gsgam2h,rgs,rgsh, ix..jz, dt)",
                "in": "state = mpF_drv_*_out; plus mpF_pt_Vt (WPREC), "
                      "mpF_pt_qd_in, mpF_pt_rgs, mpF_pt_rgsh, mpF_gsgam2, "
                      "mpF_gsgam2h, mpF_z, mpF_zh, mpF_ix..mpF_jz, dt",
                "out": "mpF_pt_*_out + mpF_pt_precip*",
                "note": "Called ONCE on the post-mp_driver state. Vt is dumped "
                        "as an input so the mini-app does not have to reproduce "
                        "mp_terminal_velocity to test the transport itself.",
            },
            "3_SATURATION_adjustment": {
                "call": "SATURATION_setrange(100.0,90.0); "
                        "SATURATION_adjustment(ijdim,kdim, rhog, rhoge, rhogq, "
                        "tem, q, qd, gsgam2, ice_adjust=.false.)",
                "in": "state = mpF_pt_*_out (rhog, rhoge, rhogq, tem, q), "
                      "mpF_pt_qd_out, mpF_gsgam2",
                "out": "mpF_sat_rhoge_out, mpF_sat_rhogq_out, mpF_sat_tem_out, "
                       "mpF_sat_q_out",
            },
        },
        "scalars": SCALARS,
        "arrays": records,
        "totals": {
            "n_arrays": len(records),
            "raw_bytes": total_bin,
        },
    }

    with open(os.path.join(size_dir, "manifest.json"), "w") as fh:
        json.dump(doc, fh, indent=2)
    if do_npy:
        with open(os.path.join(npy_dir, "manifest.json"), "w") as fh:
            json.dump(doc, fh, indent=2)

    # ---- human-readable text manifest -------------------------------------
    with open(os.path.join(size_dir, "manifest.txt"), "w") as fh:
        w = fh.write
        w(f"# NICAM NSW6 cloud-microphysics mini-app dataset -- size variant: {size}\n")
        w(f"# generated {generated} from {dump_dir}\n#\n")
        w("# PROVENANCE\n")
        for k, v in PROVENANCE.items():
            w(f"#   {k:18s} : {v}\n")
        w("#\n# DIMENSIONS\n")
        for k, v in DIMS.items():
            w(f"#   {k:6s} : {v}\n")
        w("#\n# FORMATS\n")
        for k, v in doc["formats"].items():
            w(f"#   {k:4s} : {v}\n")
        w(f"#   region layout : {doc['region_layout']}\n")
        w("#\n# ENTRY POINTS\n")
        for k, v in doc["entry_points"].items():
            w(f"#   [{k}]\n")
            for kk, vv in v.items():
                w(f"#      {kk:5s}: {vv}\n")
        w("#\n# SCALAR PARAMETERS (not stored as .bin arrays; constant for this run)\n")
        for group, gv in SCALARS.items():
            w(f"#   --- {group} ---\n")
            for k, v in gv.items():
                if k == "_note":
                    w(f"#       note: {v}\n")
                elif isinstance(v, dict):
                    val = v.get("value")
                    un = v.get("units", "")
                    mn = v.get("meaning", "")
                    w(f"#       {k:22s} = {val!r:26s} [{un}] {mn}\n")
                else:
                    w(f"#       {k:22s} = {v!r}\n")
        w("#\n")
        w("# ARRAY INVENTORY\n")
        w("# columns: name | role | shape(dims) | dtype | byte-order | units | "
          "min | max | mean | nonzero-frac | meaning\n")
        w("#\n")
        hdr = (f"{'NAME':24s} {'ROLE':14s} {'SHAPE':22s} {'DTYPE':8s} {'ENDIAN':10s} "
               f"{'UNITS':9s} {'MIN':>13s} {'MAX':>13s} {'MEAN':>13s} {'NZFRAC':>8s}  MEANING\n")
        w("# " + hdr)
        for r in records:
            shape = "(" + ",".join(f"{n}={d}" for n, d in
                                   zip(r["dim_names"], r["shape"])) + ")"
            s = r["stats"]
            w(f"{r['name']:24s} {r['role']:14s} {shape:22s} {r['dtype']:8s} "
              f"{'big':10s} {r['units']:9s} "
              f"{s['min']:13.5e} {s['max']:13.5e} {s['mean']:13.5e} "
              f"{s['nonzero_fraction']:8.4f}  {r['description']}\n")
        w("\n\n# PER-TRACER / PER-PHASE BREAKDOWN (arrays with a v or phase axis)\n")
        w(f"# {'ARRAY':24s} {'SLICE':8s} {'MIN':>13s} {'MAX':>13s} {'MEAN':>13s} {'NZFRAC':>8s}\n")
        for r in records:
            if not r.get("slice_stats"):
                continue
            for lab, s in r["slice_stats"].items():
                w(f"  {r['name']:24s} {lab:8s} {s['min']:13.5e} {s['max']:13.5e} "
                  f"{s['mean']:13.5e} {s['nonzero_fraction']:8.4f}\n")
        w(f"\n# {len(records)} arrays, {total_bin} bytes of raw records "
          f"({total_bin/2**30:.2f} GiB)\n")

    return total_bin


# --------------------------------------------------------------------------
def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # nicam_mp/
    root = os.path.dirname(here)                                        # workspace
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dump", default=os.path.join(root, "run/mp-dataset-realstep/dump"),
                    help="kksref Stage-F dump directory")
    ap.add_argument("--out", default=os.path.join(here, "dataset"),
                    help="output dataset root")
    ap.add_argument("--sizes", default="full,region0",
                    help="comma-separated size variants to build")
    ap.add_argument("--lall", type=int, default=10, help="ADM_lall of the dump")
    ap.add_argument("--no-raw", action="store_true", help="skip the .bin package")
    ap.add_argument("--no-npy", action="store_true", help="skip the .npy mirror")
    ap.add_argument("--exclude", default="",
                    help="comma-separated name PREFIXES to drop from the package "
                         "(e.g. 'mpF_pt_,mpF_sat_' to keep only the whole-step "
                         "mp_driver arrays). Applies to every size built in this "
                         "invocation, so build trimmed and untrimmed sizes in "
                         "separate invocations.")
    ap.add_argument("--manifest-only", action="store_true",
                    help="rewrite manifest.txt / manifest.json only "
                         "(no data files are written)")
    args = ap.parse_args()
    if args.manifest_only:
        args.no_raw = args.no_npy = True

    for size in [s.strip() for s in args.sizes.split(",") if s.strip()]:
        if size not in ("full", "region0"):
            sys.exit(f"ERROR: unknown size '{size}' (use full and/or region0)")
        t0 = time.time()
        exclude = tuple(p.strip() for p in args.exclude.split(",") if p.strip())
        size_dir, npy_dir, records = build(args.dump, args.out, size, args.lall,
                                           not args.no_raw, not args.no_npy,
                                           exclude)
        total = write_manifests(size_dir, npy_dir, records, size, args.dump,
                                (not args.no_npy) or args.manifest_only)
        print(f"[{size}] {len(records)} arrays, raw {total/2**30:.2f} GiB "
              f"-> {size_dir}  ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
