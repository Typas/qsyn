#!/usr/bin/env python3
"""2x2 FastTODD benchmark table: {qsyn, Vandaele} x {pyzx, raw}.

Reads bench2x2/all.csv (long format: circuit,run,n,T,seconds,status,note), or
assembles it live from bench2x2/results/*.row, and pivots on `run` into a 2x2
comparison table. Every cell reflects an actually-observed run outcome:

  <T>/<s>s  finished      TLE  timed out at cap       ERR  errored
  MLE       memory cap    N/A  input rejected (Vandaele: CZ/Tof) or unavailable
  —         not attempted (no row recorded)
"""
import csv
import glob
from pathlib import Path

ROOT = Path("/home/liao/qsyn")
PAPER = ROOT / "benchmark" / "paper-todd.csv"
B = ROOT / "bench2x2"
RUNS = ["qsyn_pyzx", "vand_pyzx", "qsyn_raw", "vand_raw"]
HEAD = {"qsyn_pyzx": "qsyn·pyzx", "vand_pyzx": "Vandaele·pyzx",
        "qsyn_raw": "qsyn·raw", "vand_raw": "Vandaele·raw"}


def fmt(t, s, status):
    st = (status or "").strip()
    if st == "ok":
        return f"{t}/{s}s" if t not in ("", "NA", None) else "ERR"
    if st in ("TLE", "MLE", "ERR"):
        return st
    if st in ("NA", "N/A"):
        return "N/A"
    return "—"


def load_rows():
    rows = []
    allcsv = B / "all.csv"
    if allcsv.exists():
        rows = list(csv.DictReader(open(allcsv)))
    else:  # live: assemble from per-job rows
        cols = ["circuit", "run", "n", "T", "seconds", "status", "note"]
        for f in sorted(glob.glob(str(B / "results" / "*.row"))):
            for rec in csv.reader(open(f)):
                if len(rec) >= 6:
                    rows.append(dict(zip(cols, rec + [""] * (len(cols) - len(rec)))))
    return rows


def main():
    rows = load_rows()
    cell, nq = {}, {}
    for r in rows:
        c, run = r["circuit"], r["run"]
        cell[(c, run)] = fmt(r.get("T"), r.get("seconds"), r.get("status"))
        if r.get("n") not in ("", "?", None) and c not in nq:
            nq[c] = r["n"]

    order = [r["circuit"] for r in csv.DictReader(open(PAPER))]
    seen = set()
    out = [
        "# 2×2 FastTODD benchmark — T-count / wall seconds",
        "",
        "Legend: `T/s`=finished · `TLE`=timed out (1800s cap) · `MLE`=memory cap · "
        "`ERR`=errored · `N/A`=input rejected (Vandaele: CZ/Tof) or unavailable · "
        "`—`=not attempted",
        "",
        "| circuit | n | " + " | ".join(HEAD[r] for r in RUNS) + " |",
        "|---|--:|" + "|".join(["---"] * len(RUNS)) + "|",
    ]
    circuits = order + [c for (c, _) in cell if c not in order]
    for c in circuits:
        if c in seen:
            continue
        seen.add(c)
        if not any((c, run) in cell for run in RUNS):
            continue
        cells = " | ".join(cell.get((c, run), "—") for run in RUNS)
        out.append(f"| {c} | {nq.get(c, '?')} | {cells} |")
    text = "\n".join(out) + "\n"
    (B / "table.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
