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
"""

import argparse
import collections
import re
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
PRIMER_RE = re.compile(r"^(?P<scheme>.+)_(?P<amp>\d+)_(?P<dir>LEFT|RIGHT)(?:_(?P<alt>\S+))?$")


def revcomp(s):
    return s.translate(COMP)[::-1]


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
                dict(start=int(f[1]), end=int(f[2]), name=f[3]))
    return amps


def apply_variants(seq, path):
    """Apply a TSV of <0-based pos>\t<ref base>\t<alt base>, verifying each ref base."""
    seq = list(seq)
    applied = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                sys.exit(f"{path}:{lineno}: expected 'pos ref alt'")
            pos, ref, alt = int(parts[0]), parts[1].upper(), parts[2].upper()
            if seq[pos] != ref:
                sys.exit(f"{path}:{lineno}: reference has {seq[pos]} at {pos}, not {ref}")
            seq[pos] = alt
            applied.append((pos, ref, alt))
    return "".join(seq), applied


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
    ap.add_argument("--write-ref", help="write the (possibly mutated) reference here")
    args = ap.parse_args()

    refname, ref = load_ref(args.ref)
    if args.variants:
        ref, applied = apply_variants(ref, args.variants)
        print(f"applied {len(applied)} substitution(s): "
              + ", ".join(f"{p}{r}>{a}" for p, r, a in applied))
    if args.write_ref:
        with open(args.write_ref, "w") as fh:
            fh.write(f">{refname}\n")
            for i in range(0, len(ref), 60):
                fh.write(ref[i:i + 60] + "\n")

    amps = load_primers(args.bed)
    pairs = 0
    short = []
    with open(args.r1, "w") as r1, open(args.r2, "w") as r2:
        for amp in sorted(amps):
            left, right = amps[amp]["LEFT"], amps[amp]["RIGHT"]
            if not left or not right:
                sys.exit(f"amplicon {amp} is missing a LEFT or RIGHT primer")
            for lp in left:
                for rp in right:
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

    print(f"ref={refname} len={len(ref)} amplicons={len(amps)} pairs={pairs}")
    if short:
        print(f"WARNING: {len(sorted(set(short)))} amplicon(s) shorter than "
              f"--read-len {args.read_len}, skipped: {sorted(set(short))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
