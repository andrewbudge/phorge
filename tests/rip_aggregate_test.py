#!/usr/bin/env python3
"""Checks for `rip --aggregate` (one multi-FASTA per canonical gene)."""
import argparse
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run(c_bin: str, argv: list[str]) -> tuple[int, bytes, bytes]:
    p = subprocess.run([c_bin, *argv], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.returncode, p.stdout, p.stderr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--c", default=str(ROOT / "build" / "phorge"))
    args = ap.parse_args()

    fixture = ROOT / "tests" / "fixtures" / "seqripper" / "NC_001700.gb"
    work = ROOT / "tests" / "_out" / "rip-aggregate"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True, exist_ok=True)
    fail = 0

    def check(name, ok, detail=""):
        nonlocal fail
        print(f"[{'OK' if ok else 'FAIL':4s}] {name}" + (f"  {detail}" if detail and not ok else ""))
        if not ok:
            fail += 1

    out = work / "all"
    rc, _, se = run(args.c, ["rip", "--all", "--no-trna", "--aggregate", "-o", str(out), str(fixture)])
    names = {p.name for p in out.iterdir()}
    check("rip --aggregate exit", rc == 0, f"rc={rc} err={se[:200]!r}")
    check("rip --aggregate one file per gene", len(names) == 16 and "COI.fasta" in names,
          f"names={sorted(names)}")
    coi = (out / "COI.fasta").read_text()
    check("rip --aggregate collapses gene/CDS copies", coi.count(">") == 1,
          f"records={coi.count('>')}")
    log = (out / "seqripper_log.tsv").read_text()
    check("rip --aggregate canonical gene names", "NC_001700.1\tCOI\t" in log
          and "NC_001700.1\tCOX1\t" not in log)

    variant = work / "variant.gb"
    variant.write_text(fixture.read_text().replace("VERSION     NC_001700.1",
                                                   "VERSION     TEST00001.1")
                                      .replace('/gene="COX1"', '/gene="CO1"'))
    out_many = work / "many"
    rc, _, se = run(args.c, ["rip", "-g", "COI", "-a", "-o", str(out_many),
                             str(fixture), str(variant)])
    coi_many = (out_many / "COI.fasta").read_text()
    check("rip -a spans genomes and synonyms", rc == 0 and coi_many.count(">") == 2,
          f"rc={rc} records={coi_many.count('>')}")

    out_miss = work / "missing"
    rc, _, _ = run(args.c, ["rip", "-g", "COI,NOTAGENE", "--aggregate", "-o", str(out_miss),
                            str(fixture)])
    log_miss = (out_miss / "seqripper_log.tsv").read_text()
    check("rip --aggregate reports missing genes",
          rc == 0 and "\tNOTAGENE\t0\tmissing\t-\n" in log_miss
          and sorted(p.name for p in out_miss.iterdir()) == ["COI.fasta", "seqripper_log.tsv"])

    synth = work / "coiii.gb"
    synth.write_text('VERSION     SYNTH1.1\n'
                     'SOURCE      mitochondrion syntheticus\n'
                     '  ORGANISM  syntheticus\n'
                     'FEATURES             Location/Qualifiers\n'
                     '     CDS             1..9\n'
                     '                     /product="cytochrome c oxidase subunit III"\n'
                     'ORIGIN\n'
                     '        1 atgcatgcat\n'
                     '//\n')
    out_roman = work / "roman"
    rc, _, _ = run(args.c, ["rip", "-g", "COIII", "--aggregate", "-o", str(out_roman), str(synth)])
    check("rip --aggregate keeps subunit III distinct from I",
          rc == 0 and (out_roman / "COIII.fasta").exists()
          and (out_roman / "COIII.fasta").read_text().count(">") == 1)

    print(f"rip aggregate cases checked; {fail} failure(s)")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
