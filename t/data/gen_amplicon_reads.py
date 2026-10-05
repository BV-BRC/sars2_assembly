#!/usr/bin/env python3
"""Generate synthetic paired FASTQ reads for every amplicon in a primer BED.

Used by t/server-tests/ivar-trim-roundtrip.t to exercise the alignment and
primer-trimming path with data whose correct output is known exactly.

Reads are emitted per LEFT/RIGHT primer *combination*, not per amplicon, so
alternate primers at a shared binding site are each exercised rather than one
arbitrarily winning. Primer names are parsed by anchoring on _LEFT/_RIGHT; a
positional split on "_" does not work, because scheme names themselves contain
underscores and alt suffixes vary (_0, _1, _alt1).

    gen_amplicon_reads.py --ref R.fasta --bed P.bed --r1 R1.fq --r2 R2.fq

With --variants, substitutions are applied to the reference before reads are
cut, so a mutation can be placed under a primer-binding site to check that it
is still recovered from the overlapping neighbouring amplicon after trimming.

With --anneal-model, an amplicon is only amplified when at least one LEFT and
one RIGHT oligo actually anneals to the (possibly mutated) template, so a
primer-site mutation can knock the amplicon out entirely. This is what makes a
dropout test possible: without it every amplicon amplifies unconditionally and
a scheme that cannot bind the template looks identical to one that can.

The model is a deliberate simplification of PCR and is the one assumption in
this test that is not a measured fact -- see --anneal-model for what each rule
claims.
"""

import argparse
import collections
import re
import sys

COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                     "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
PRIMER_RE = re.compile(r"^(?P<scheme>.+)_(?P<amp>\d+)_(?P<dir>LEFT|RIGHT)(?:_(?P<alt>\S+))?$")

# Degenerate oligo bases match any of their constituent bases. V5.3.2 ships one
# such primer (_84_RIGHT_2, an R) standing in for the two explicit oligos that
# upstream now ships separately; treating R as a literal would report it as a
# permanent mismatch against every template.
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def revcomp(s):
    return s.translate(COMP)[::-1]


def mismatches(oligo, template):
    """Mismatch positions of an oligo against the template it would anneal to.

    Both are given 5'->3' in the oligo's own orientation, so index 0 is the 5'
    end and index len-1 is the 3' end -- the base the polymerase extends from.
    Returns [(dist_from_3prime, oligo_base, template_base), ...].
    """
    out = []
    for i, (o, t) in enumerate(zip(oligo, template)):
        if t not in IUPAC.get(o, o):
            out.append((len(oligo) - 1 - i, o, t))
    return out


def load_ref(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    sys.exit(f"{path}: expected a single FASTA record")
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    if name is None:
        sys.exit(f"{path}: no FASTA record found")
    return name, "".join(seq).upper()


def load_primers(path):
    amps = collections.defaultdict(lambda: {"LEFT": [], "RIGHT": []})
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 4:
                sys.exit(f"{path}:{lineno}: fewer than 4 columns")
            m = PRIMER_RE.match(f[3])
            if not m:
                sys.exit(f"{path}:{lineno}: cannot parse primer name {f[3]!r}")
            amps[int(m.group("amp"))][m.group("dir")].append(
                dict(start=int(f[1]), end=int(f[2]), name=f[3],
                     strand=f[5] if len(f) > 5 else None,
                     seq=f[6].upper() if len(f) > 6 else None))
    return amps


def anneals(primer, template, model, max_mismatch, window):
    """Decide whether an oligo primes synthesis on this template.

    Returns (ok, mismatch_list). The mismatch list is always computed, so a
    caller can report *why* an amplicon dropped out even under --anneal-model
    none.
    """
    site = template[primer["start"]:primer["end"]]
    # The BED sequence column is already written in the oligo's own 5'->3'
    # orientation: as the top strand for a LEFT primer, as its reverse
    # complement for a RIGHT primer. So the template must be flipped to match
    # for RIGHT, and in both cases index len-1 is the 3' end.
    mm = mismatches(primer["seq"], site if primer["strand"] != "-" else revcomp(site))
    if model == "none":
        return True, mm
    if model == "strict":
        return not mm, mm
    # three-prime: a mismatch in the last `window` bases blocks extension
    # outright, and two or more mismatches anywhere depress Tm enough to lose
    # the amplicon. A single mismatch further 5' is tolerated -- in the lab it
    # reduces yield rather than abolishing it, which this binary model cannot
    # express. See the --anneal-model help text.
    if any(d3 < window for d3, _, _ in mm):
        return False, mm
    return len(mm) <= max_mismatch, mm


def load_scenario(path, name):
    """Read the named row set out of the 4-column scenario table."""
    subs, seen = [], set()
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                sys.exit(f"{path}:{lineno}: expected 'set pos ref alt'")
            seen.add(parts[0])
            if parts[0] == name or name == "all":
                subs.append((int(parts[1]), parts[2].upper(), parts[3].upper()))
    if not subs and name != "none":
        sys.exit(f"{path}: no scenario named {name!r} (have: {', '.join(sorted(seen))})")
    return subs


def load_variants(path):
    """Read a TSV of <0-based pos>\t<ref base>\t<alt base>."""
    subs = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                sys.exit(f"{path}:{lineno}: expected 'pos ref alt'")
            subs.append((int(parts[0]), parts[1].upper(), parts[2].upper()))
    return subs


def apply_variants(seq, subs):
    """Apply substitutions, verifying the stated reference base at each site."""
    seq = list(seq)
    for pos, ref, alt in subs:
        if seq[pos] != ref:
            sys.exit(f"reference has {seq[pos]} at 0-based {pos}, not {ref}")
        seq[pos] = alt
    return "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="reference FASTA (single record)")
    ap.add_argument("--bed", required=True, help="primer BED")
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--read-len", type=int, default=250,
                    help="read length; must be long enough to tile the insert (default 250)")
    ap.add_argument("--depth", type=int, default=6,
                    help="read pairs per LEFT/RIGHT combination (default 6)")
    ap.add_argument("--variants",
                    help="TSV of '<0-based pos> <ref> <alt>' to apply before cutting reads")
    ap.add_argument("--scenarios", help="4-column scenario table: 'set pos ref alt'")
    ap.add_argument("--scenario", default="none",
                    help="row set to select from --scenarios ('all' for every row)")
    ap.add_argument("--anneal-model", choices=("none", "strict", "three-prime"),
                    default="none",
                    help="none: every amplicon amplifies regardless of primer/template "
                         "mismatch (the historical behaviour). strict: any mismatch "
                         "kills the oligo. three-prime (recommended): a mismatch within "
                         "--three-prime-window of the 3' end blocks extension, and more "
                         "than --max-mismatch mismatches anywhere kills the oligo. "
                         "An amplicon amplifies if ANY of its LEFT and ANY of its RIGHT "
                         "oligos anneal, which is what makes a spiked-in alternate "
                         "primer able to rescue a dropout.")
    ap.add_argument("--max-mismatch", type=int, default=1,
                    help="mismatches tolerated outside the 3' window (default 1)")
    ap.add_argument("--three-prime-window", type=int, default=5,
                    help="bases from the 3' end where any mismatch is fatal (default 5)")
    ap.add_argument("--dropout-report",
                    help="write a TSV of every amplicon and whether it amplified")
    ap.add_argument("--write-ref", help="write the (possibly mutated) reference here")
    args = ap.parse_args()

    refname, ref = load_ref(args.ref)
    subs = []
    if args.variants:
        subs += load_variants(args.variants)
    if args.scenarios:
        subs += load_scenario(args.scenarios, args.scenario)
    if subs:
        ref = apply_variants(ref, subs)
        print(f"applied {len(subs)} substitution(s): "
              + ", ".join(f"{p}{r}>{a}" for p, r, a in subs))
    if args.write_ref:
        with open(args.write_ref, "w") as fh:
            fh.write(f">{refname}\n")
            for i in range(0, len(ref), 60):
                fh.write(ref[i:i + 60] + "\n")

    amps = load_primers(args.bed)
    if args.anneal_model != "none":
        missing = [p["name"] for v in amps.values() for d in ("LEFT", "RIGHT")
                   for p in v[d] if not p["seq"]]
        if missing:
            sys.exit(f"--anneal-model needs the BED sequence column; "
                     f"{len(missing)} primer(s) have none, e.g. {missing[0]}")

    pairs = 0
    short, dropped, report = [], [], []
    with open(args.r1, "w") as r1, open(args.r2, "w") as r2:
        for amp in sorted(amps):
            left, right = amps[amp]["LEFT"], amps[amp]["RIGHT"]
            if not left or not right:
                sys.exit(f"amplicon {amp} is missing a LEFT or RIGHT primer")

            ok = {}
            for d, ps in (("LEFT", left), ("RIGHT", right)):
                for p in ps:
                    ok[p["name"]] = anneals(p, ref, args.anneal_model,
                                            args.max_mismatch, args.three_prime_window)
            good_l = [p for p in left if ok[p["name"]][0]]
            good_r = [p for p in right if ok[p["name"]][0]]

            why = "; ".join(
                f"{p['name']}:" + (",".join(f"{o}>{t}@{d3}nt-from-3'"
                                            for d3, o, t in ok[p["name"]][1]) or "match")
                for p in left + right if ok[p["name"]][1])
            report.append((amp, "AMPLIFIED" if (good_l and good_r) else "DROPOUT",
                           f"{len(good_l)}/{len(left)}", f"{len(good_r)}/{len(right)}",
                           why or "-"))
            if not (good_l and good_r):
                dropped.append(amp)
                continue

            for lp in good_l:
                for rp in good_r:
                    seq = ref[lp["start"]:rp["end"]]
                    if len(seq) < args.read_len:
                        short.append(amp)
                        continue
                    fwd = seq[:args.read_len]
                    rev = revcomp(seq[-args.read_len:])
                    for _ in range(args.depth):
                        pairs += 1
                        tag = f"amp{amp}_{lp['name']}__{rp['name']}_{pairs}"
                        r1.write(f"@{tag}/1\n{fwd}\n+\n{'I' * len(fwd)}\n")
                        r2.write(f"@{tag}/2\n{rev}\n+\n{'I' * len(rev)}\n")

    if args.dropout_report:
        with open(args.dropout_report, "w") as fh:
            fh.write("amplicon\tstatus\tleft_ok\tright_ok\tmismatches\n")
            for row in report:
                fh.write("\t".join(str(c) for c in row) + "\n")

    print(f"ref={refname} len={len(ref)} amplicons={len(amps)} pairs={pairs} "
          f"model={args.anneal_model}")
    if dropped:
        print(f"DROPOUT: {len(dropped)} amplicon(s) did not amplify: {dropped}")
    if short:
        print(f"WARNING: {len(sorted(set(short)))} amplicon(s) shorter than "
              f"--read-len {args.read_len}, skipped: {sorted(set(short))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
