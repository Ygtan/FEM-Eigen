#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CTest entry point for the 3D continuous-box-beam benchmark.

Runs the current build's ``example_3d_report.exe`` on each archived .dat,
parses <stem>.report.txt, and compares against frozen ANSYS BEAM4 values
(Euler-Bernoulli, same theory as FEM-Eigen 3D). BEAM188 is Timoshenko and
is archived/quantified separately, not used as a golden reference.

All runtime output is ASCII so CTest logs, PowerShell and redirected files
stay readable regardless of host console encoding.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANSYS_DEFAULT = Path(r"C:\Program Files\ANSYS Inc\v180\ANSYS\bin\winx64\ANSYS180.exe")
MSYS_BIN = r"C:\msys64\mingw64\bin"

ARCHIVE = {
    "dat": "{stem}.dat",
    "inp": "{stem}.inp",
    "inp_b4": "{stem}_b4.inp",
    "model": "{stem}.model.json",
    "report": "{stem}.report.txt",
    "b4_disp": "{stem}_b4_disp.csv",
    "b4_reac": "{stem}_b4_reac.csv",
    "b4_elem": "{stem}_b4_elem.csv",
    "b188_disp": "{stem}_disp.csv",
    "b188_reac": "{stem}_reac.csv",
    "b188_elem": "{stem}_elem.csv",
    "b4_out": "{stem}_b4_ansys.out",
    "b188_out": "{stem}_ansys.out",
}
FORCE_KEYS = ("N", "Vy", "Vz", "T", "My", "Mz")


def read_table(path: Path) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        t = raw.replace(",", " ").split()
        if len(t) < 2:
            continue
        try:
            v = [float(x) for x in t]
        except ValueError:
            continue
        out[int(round(v[0]))] = v[1:]
    return out


def parse_report(path: Path) -> dict:
    """Parse example_3d_report.exe <stem>.report.txt (3D, 6 DOF/node)."""
    disp: dict[int, list[float]] = {}
    force: dict[int, dict[str, float]] = {}
    section = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "节点位移" in line:
            section = "disp"
            continue
        if "单元内力" in line:
            section = "force"
            continue
        if line.startswith("单元应力") or line.startswith("支座反力"):
            section = ""
            continue
        t = line.replace(",", " ").split()
        try:
            key = int(t[0])
        except (ValueError, IndexError):
            continue
        if section == "disp":
            try:
                v = [float(x) for x in t[1:7]]
            except ValueError:
                continue
            if len(v) == 6:
                disp[key] = v
        elif section == "force":
            kv = {}
            for tok in t[2:]:
                if "=" in tok:
                    k, _, v = tok.partition("=")
                    try:
                        kv[k] = float(v)
                    except ValueError:
                        pass
            if any(k in kv for k in FORCE_KEYS):
                force[key] = {k: kv.get(k, 0.0) for k in FORCE_KEYS}
    return {"disp": disp, "force": force}


def archive_checks(cases: Path, stem: str, doc: dict) -> list[str]:
    errs: list[str] = []
    for key, pat in ARCHIVE.items():
        p = cases / pat.format(stem=stem)
        if not p.exists():
            errs.append("missing %s (%s)" % (p.name, key))
    model_p = cases / ARCHIVE["model"].format(stem=stem)
    if model_p.exists():
        try:
            model = json.loads(model_p.read_text(encoding="utf-8"))
            nn, ne = len(model["nodes"]), len(model["elements"])
            if nn != doc.get("n_nodes") or ne != doc.get("n_elements"):
                errs.append("model size %d/%d != reference %s/%s" %
                            (nn, ne, doc.get("n_nodes"), doc.get("n_elements")))
            for key, want in (("b4_disp", nn), ("b4_reac", nn), ("b4_elem", ne),
                              ("b188_disp", nn), ("b188_reac", nn),
                              ("b188_elem", ne)):
                p = cases / ARCHIVE[key].format(stem=stem)
                if p.exists():
                    got = len(read_table(p))
                    if got != want:
                        errs.append("%s has %d rows, expect %d" %
                                    (p.name, got, want))
        except (ValueError, KeyError) as exc:
            errs.append("bad model json: %s" % exc)
    return errs


def run_fem(bin_dir: Path, cases: Path, stem: str, tmp: Path) -> tuple[dict | None, list[str]]:
    exe = bin_dir / "example_3d_report.exe"
    if not exe.exists():
        return None, ["missing solver %s" % exe]
    dat = cases / (stem + ".dat")
    work = tmp / stem
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy2(dat, work / dat.name)
    env = os.environ.copy()
    env["PATH"] = MSYS_BIN + os.pathsep + env.get("PATH", "")
    try:
        p = subprocess.run([str(exe), dat.name], cwd=str(work), env=env,
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=300)
    except Exception as exc:                                   # noqa: BLE001
        return None, ["solver exception: %s" % exc]
    if p.returncode != 0:
        return None, ["solver exit=%d: %s" % (p.returncode, (p.stderr or "")[:200])]
    rep = work / (stem + ".report.txt")
    if not rep.exists():
        return None, ["solver produced no report.txt"]
    return parse_report(rep), []


def extract(parsed: dict, chk: dict, b4reac: dict) -> float | None:
    kind = chk["kind"]
    try:
        if kind == "disp":
            return parsed["disp"][chk["node"]][chk["dof"] - 1]
        if kind == "force":
            v = parsed["force"][chk["elem"]][chk["col"]]
            return abs(v) if chk.get("abs") else v
        if kind == "reaction":
            return b4reac[chk["node"]][chk["dof"] - 1]
    except (KeyError, IndexError):
        return None
    return None


def compare_one(bin_dir: Path, cases: Path, stem: str, doc: dict,
                tmp: Path) -> tuple[list[str], list[str], int, int]:
    errs = archive_checks(cases, stem, doc)
    lines: list[str] = []
    if errs:
        return errs, lines, 0, 0
    parsed, rerr = run_fem(bin_dir, cases, stem, tmp)
    errs.extend(rerr)
    if parsed is None:
        return errs, lines, 0, 0
    b4reac = read_table(cases / ARCHIVE["b4_reac"].format(stem=stem))

    npass = nfail = 0
    for chk in doc["checks"]:
        got = extract(parsed, chk, b4reac)
        ref = float(chk["ref"])
        name = "%s:%s" % (stem, chk["name"])
        if got is None:
            lines.append("[FAIL] %-38s missing value" % name)
            nfail += 1
            continue
        if chk.get("abs"):
            got = abs(got)
        d = abs(got - ref)
        tol = float(chk.get("atol", 1e-9)) + float(chk.get("rtol", 1e-3)) * abs(ref)
        if d <= tol:
            # Near-zero reference: relative error is meaningless, show abs instead
            tag = ("rel=%-9.3g" % (d / abs(ref))) if abs(ref) > 1e-6 \
                else ("abs=%-9.3g" % d)
            lines.append("[PASS] %-38s got=%-15.8g ref=%-15.8g %s" %
                         (name, got, ref, tag))
            npass += 1
        else:
            lines.append("[FAIL] %-38s got=%-15.8g ref=%-15.8g diff=%-9.3g tol=%g" %
                         (name, got, ref, d, tol))
            nfail += 1
    if nfail:
        errs.append("%s has %d failed checks" % (stem, nfail))
    return errs, lines, npass, nfail


def live_ansys(cases: Path, stem: str, ansys: Path,
               tmp: Path) -> tuple[list[str], list[str]]:
    """Optional: re-run both APDL jobs and compare with frozen ANSYS CSV."""
    errs: list[str] = []
    lines: list[str] = []
    if not ansys.exists():
        return ["ANSYS not found: %s" % ansys], lines
    work = tmp / (stem + "_live")
    work.mkdir(parents=True, exist_ok=True)
    # Each job is compared against its OWN frozen CSV (BEAM188 vs BEAM188,
    # BEAM4 vs BEAM4); cross-theory difference is documented separately.
    jobs = [
        ("%s.inp" % stem, stem,
         [("%s_disp.csv" % stem, "%s_disp.csv" % stem),
          ("%s_reac.csv" % stem, "%s_reac.csv" % stem),
          ("%s_elem.csv" % stem, "%s_elem.csv" % stem)]),
        ("%s_b4.inp" % stem, stem + "_b4",
         [("%s_b4_disp.csv" % stem, "%s_b4_disp.csv" % stem),
          ("%s_b4_reac.csv" % stem, "%s_b4_reac.csv" % stem),
          ("%s_b4_elem.csv" % stem, "%s_b4_elem.csv" % stem)]),
    ]
    for inp_name, job, pairs in jobs:
        shutil.copy2(cases / inp_name, work / inp_name)
        try:
            p = subprocess.run([str(ansys), "-b", "-j", job, "-i", inp_name,
                                "-o", job + ".out"], cwd=str(work),
                               capture_output=True, text=True, encoding="utf-8",
                               errors="replace", timeout=600)
        except Exception as exc:                               # noqa: BLE001
            errs.append("ANSYS %s exception: %s" % (job, exc))
            continue
        if p.returncode != 0:
            errs.append("ANSYS %s exit=%d" % (job, p.returncode))
            continue
        for produced, frozen in pairs:
            pf, pz = work / produced, cases / frozen
            if not pf.exists() or not pz.exists():
                errs.append("LIVE %s missing %s" % (job, produced))
                continue
            a, b = read_table(pf), read_table(pz)
            if set(a) != set(b):
                errs.append("LIVE %s row keys differ" % produced)
                continue
            dmax = 0.0
            scale = 0.0
            for k in a:
                for x, y in zip(a[k], b[k]):
                    dmax = max(dmax, abs(x - y))
                    scale = max(scale, abs(y))
            rel = dmax / (scale or 1.0)
            ok = rel <= 1e-3
            lines.append("  LIVE %-34s maxdiff=%-10.3g rel=%-9.3g %s" %
                         (frozen, dmax, rel, "PASS" if ok else "FAIL"))
            if not ok:
                errs.append("LIVE %s out of tolerance" % frozen)
    return errs, lines


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--cases", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--live-ansys", action="store_true")
    ap.add_argument("--ansys", type=Path, default=ANSYS_DEFAULT)
    a = ap.parse_args()

    docs = []
    for p in sorted(a.reference.glob("*_ansys.json")):
        docs.append(json.loads(p.read_text(encoding="utf-8")))
    if not docs:
        print("[FAIL] no reference JSON in %s" % a.reference)
        return 1

    all_errs: list[str] = []
    total_pass = total_fail = 0
    with tempfile.TemporaryDirectory(prefix="cb3d_") as td:
        tmp = Path(td)
        for doc in docs:
            stem = doc["id"]
            print("== %s (%s, %d checks) ==" %
                  (stem, doc["dimension"], len(doc["checks"])))
            errs, lines, np_, nf_ = compare_one(a.bin, a.cases, stem, doc, tmp)
            print("\n".join(lines))
            total_pass += np_
            total_fail += nf_
            if errs:
                for e in errs:
                    print("[FAIL] %s" % e)
                all_errs.extend(errs)
            if a.live_ansys:
                print("-- LIVE ANSYS cross-check: %s --" % stem)
                lerr, llines = live_ansys(a.cases, stem, a.ansys, tmp)
                print("\n".join(llines))
                for e in lerr:
                    print("[FAIL] %s" % e)
                all_errs.extend(lerr)
    print("Result: %d/%d checks passed, %d archive/live error(s)" %
          (total_pass, total_pass + total_fail, len(all_errs)))
    return 1 if (all_errs or total_fail) else 0


if __name__ == "__main__":
    sys.exit(main())
