"""
Deep dive into the most promising Level 6 candidates.

Initial exploration revealed:
  - Third base T and C give IDENTICAL side chain sums (864 each!)
  - Third base M=(T,G) side chains = 1776 = SM(999+777)
  - Third base K=(C,A) full weight(stop=74) = 3996 = SM(3×999+999)
  - GC=1 codons: side chains = 1332 = SM(999+333)
  - YYR - RRY pattern difference = 74 (the constant part!)
  - All reverse complement pair sums = 3404 = 37×92 = 74×46
"""

from math import gcd
from collections import defaultdict
from itertools import combinations

GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'stop', 'TAG': 'stop',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'stop', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AMINO_ACIDS = {
    'G': 1, 'A': 15, 'S': 31, 'P': 42, 'V': 43, 'T': 45, 'C': 47,
    'L': 57, 'I': 57, 'N': 58, 'D': 59, 'Q': 72, 'K': 72,
    'E': 73, 'M': 75, 'H': 81, 'F': 91, 'R': 100, 'Y': 107, 'W': 130,
}

FULL_WEIGHT = {
    'G': 75, 'A': 89, 'S': 105, 'P': 115, 'V': 117, 'T': 119, 'C': 121,
    'L': 131, 'I': 131, 'N': 132, 'D': 133, 'Q': 146, 'K': 146,
    'E': 147, 'M': 149, 'H': 155, 'F': 165, 'R': 174, 'Y': 181, 'W': 204,
}

BASES = ['T', 'C', 'A', 'G']


def sc(aa):
    if aa == 'stop': return 0
    return 41 if aa == 'P' else AMINO_ACIDS[aa]


def fw(aa, stop=0):
    if aa == 'stop': return stop
    return FULL_WEIGHT[aa]


def is_sm(n):
    if n <= 0: return False
    for m in range(n // 999 + 1):
        rem = n - m * 999
        if rem >= 0 and rem % 111 == 0 and rem // 111 <= 9:
            return True
    return False


def sm_repr(n):
    if n <= 0: return str(n)
    for m in range(n // 999 + 1):
        rem = n - m * 999
        if rem >= 0 and rem % 111 == 0:
            k = rem // 111
            if 0 <= k <= 9:
                if m == 0: return f"{k}×111"
                return f"{m}×999+{k}×111" if k > 0 else f"{m}×999"
    return str(n)


def props(n):
    """One-line property summary."""
    parts = []
    if is_sm(n): parts.append(f"SM({sm_repr(n)})")
    if n > 0 and n % 37 == 0: parts.append(f"37×{n // 37}")
    if n > 0 and n % 74 == 0: parts.append(f"74×{n // 74}")
    if n > 0 and n % 111 == 0: parts.append(f"111×{n // 111}")
    return ' '.join(parts) if parts else ''


def sm_products(columns):
    products = []
    for b1, b2 in columns:
        seen = set()
        for b3 in BASES:
            aa = GENETIC_CODE[b1 + b2 + b3]
            if aa not in seen:
                seen.add(aa)
                products.append(aa)
    return products


def is_homogeneous(b1, b2):
    return len(set(GENETIC_CODE[b1 + b2 + b3] for b3 in BASES)) == 1


# =============================================================================

def finding_1_third_base_symmetry():
    """FINDING: Third base T and C produce IDENTICAL side chain sums."""
    print("\n" + "=" * 70)
    print("FINDING 1: THIRD BASE WOBBLE SYMMETRY")
    print("=" * 70)

    for b3 in BASES:
        codons = [c for c in GENETIC_CODE if c[2] == b3]
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        total_fw0 = sum(fw(GENETIC_CODE[c], stop=0) for c in codons)
        total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        print(f"  Third base = {b3}: sc={total_sc}  fw(stop=0)={total_fw0}  fw(stop=74)={total_fw74}")

    print()
    # The key finding: T and C give identical sums!
    t_sc = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE if c[2] == 'T')
    c_sc = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE if c[2] == 'C')
    print(f"  ★ Third base T side chains = {t_sc}")
    print(f"  ★ Third base C side chains = {c_sc}")
    print(f"  ★ EQUAL: {t_sc} = {c_sc} = {t_sc == c_sc}!")
    print(f"  ★ 864 = {864 // 32}×32 = {864 // 16}×16 = {864 // 12}×12")
    print(f"  ★ 864 / 37 = {864 / 37:.4f} (not divisible)")

    # Is this equality forced by the structure of the code?
    # The third base Y (T or C) distinguishes codons within 2-fold degenerate families.
    # For homogeneous columns (4-fold degenerate), T and C encode the same product.
    # For heterogeneous columns, the Y/R split determines the product.
    # The equality means: sum over all columns of [product encoded by ...T] = [product encoded by ...C]
    print("\n  Column-by-column breakdown:")
    total_t, total_c = 0, 0
    for b1 in BASES:
        for b2 in BASES:
            aa_t = GENETIC_CODE[b1 + b2 + 'T']
            aa_c = GENETIC_CODE[b1 + b2 + 'C']
            st, scc = sc(aa_t), sc(aa_c)
            total_t += st
            total_c += scc
            marker = "" if aa_t == aa_c else f" ← DIFFER (Δ={st - scc:+d})"
            print(f"    {b1}{b2}: ...T→{aa_t}({st})  ...C→{aa_c}({scc}){marker}")
    print(f"  Total: T={total_t}, C={total_c}")

    # Count differing columns
    diffs = []
    for b1 in BASES:
        for b2 in BASES:
            aa_t = GENETIC_CODE[b1 + b2 + 'T']
            aa_c = GENETIC_CODE[b1 + b2 + 'C']
            if aa_t != aa_c:
                diffs.append((b1 + b2, aa_t, aa_c, sc(aa_t) - sc(aa_c)))

    print(f"\n  Columns where T≠C: {len(diffs)}")
    for col, t, c, d in diffs:
        print(f"    {col}: {t}({sc(GENETIC_CODE[col + 'T'])}) vs {c}({sc(GENETIC_CODE[col + 'C'])}) Δ={d:+d}")
    total_delta = sum(d for _, _, _, d in diffs)
    print(f"  Sum of deltas: {total_delta} (must be 0 for equality)")

    # Now do the same for A vs G
    a_sc = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE if c[2] == 'A')
    g_sc = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE if c[2] == 'G')
    print(f"\n  Third base A side chains = {a_sc}")
    print(f"  Third base G side chains = {g_sc}")
    print(f"  Difference A-G: {a_sc - g_sc}")

    # Also check full weights
    t_fw = sum(fw(GENETIC_CODE[c], stop=0) for c in GENETIC_CODE if c[2] == 'T')
    c_fw = sum(fw(GENETIC_CODE[c], stop=0) for c in GENETIC_CODE if c[2] == 'C')
    print(f"\n  Full weight (stop=0): T={t_fw}, C={c_fw}, EQUAL={t_fw == c_fw}")
    t_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in GENETIC_CODE if c[2] == 'T')
    c_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in GENETIC_CODE if c[2] == 'C')
    print(f"  Full weight (stop=74): T={t_fw74}, C={c_fw74}, EQUAL={t_fw74 == c_fw74}")


def finding_2_third_base_rumer_sm():
    """FINDING: Third base Rumer pairs give SM numbers."""
    print("\n" + "=" * 70)
    print("FINDING 2: THIRD BASE RUMER PAIRS → SM NUMBERS")
    print("=" * 70)

    for label, bases in [("M=(T,G)", "TG"), ("K=(C,A)", "CA")]:
        codons = [c for c in GENETIC_CODE if c[2] in bases]
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        total_fw0 = sum(fw(GENETIC_CODE[c], stop=0) for c in codons)
        total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        print(f"\n  Third base ∈ {label} (32 codons):")
        print(f"    Side chains:       {total_sc}  {props(total_sc)}")
        print(f"    Full weight (s=0): {total_fw0}  {props(total_fw0)}")
        print(f"    Full weight (s=74):{total_fw74}  {props(total_fw74)}")

    # The paper found: first base Y/R and M/K patterns. Third base M/K is NEW.
    # Third base M: sc=1776 = 999+777 (SM!) and fw74=4144 = 37×112 = 74×56
    # Third base K: sc=1628 = 37×44 = 74×22 and fw74=3996 = 3×999+999 (SM!)
    print("\n  ★ Third base M: sc = 1776 = 999+777  (SM number!)")
    print("  ★ Third base K: fw(stop=74) = 3996 = 3×999+999  (SM number!)")
    print("  ★ Both third-base Rumer groups have ×37 and ×74 divisibility")
    print(f"  ★ sc_M + sc_K = {1776 + 1628} = {props(1776 + 1628)}")
    print(f"  ★ fw74_M + fw74_K = {4144 + 3996} = {props(4144 + 3996)}")
    print(f"  ★ GCD(sc_M, sc_K) = {gcd(1776, 1628)} = {props(gcd(1776, 1628))}")
    print(f"  ★ GCD(fw74_M, fw74_K) = {gcd(4144, 3996)} = {props(gcd(4144, 3996))}")


def finding_3_gc_content():
    """FINDING: GC=1 codons are deeply SM."""
    print("\n" + "=" * 70)
    print("FINDING 3: GC-CONTENT GROUPINGS")
    print("=" * 70)

    gc_groups = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        gc = sum(1 for b in codon if b in 'GC')
        gc_groups[gc].append((codon, aa))

    for gc in sorted(gc_groups):
        codons = gc_groups[gc]
        total_sc = sum(sc(aa) for _, aa in codons)
        total_fw0 = sum(fw(aa) for _, aa in codons)
        total_fw74 = sum(fw(aa, stop=74) for _, aa in codons)
        n = len(codons)
        print(f"\n  GC={gc} ({n} codons):")
        print(f"    Side chains:        {total_sc}  {props(total_sc)}")
        print(f"    Full weight (s=0):  {total_fw0}  {props(total_fw0)}")
        print(f"    Full weight (s=74): {total_fw74}  {props(total_fw74)}")

    # Key: GC=1 is spectacularly SM
    print("\n  ★ GC=1 (24 codons): sc = 1332 = 999+333 (SM!)")
    print("  ★ GC=1 (24 codons): fw = 3108 = 3×999+111 (SM!)")
    print("  ★ GC=0 + GC=1 = purely AT or 1 GC:")
    at_heavy_sc = sum(sc(aa) for gc in [0, 1] for _, aa in gc_groups[gc])
    gc_heavy_sc = sum(sc(aa) for gc in [2, 3] for _, aa in gc_groups[gc])
    print(f"    AT-heavy sc = {at_heavy_sc}  {props(at_heavy_sc)}")
    print(f"    GC-heavy sc = {gc_heavy_sc}  {props(gc_heavy_sc)}")
    print(f"    Total = {at_heavy_sc + gc_heavy_sc}  {props(at_heavy_sc + gc_heavy_sc)}")

    # What about constant parts?
    gc1_cp = sum(74 if aa != 'stop' else 0 for _, aa in gc_groups[1])
    gc1_cp74 = sum(74 for _, aa in gc_groups[1])
    print(f"\n  GC=1 constant parts (stop=0): {gc1_cp} {props(gc1_cp)}")
    print(f"  GC=1 constant parts (stop=74): {gc1_cp74} {props(gc1_cp74)}")
    # GC=1 has 24 codons, 1 is a stop (TGA has GC=1? T=0,G=1,A=0 → yes!)
    gc1_stops = [c for c, aa in gc_groups[1] if aa == 'stop']
    print(f"  GC=1 stop codons: {gc1_stops}")


def finding_4_ry_pattern_74():
    """FINDING: YYR - RRY = 74 (the constant part weight!)"""
    print("\n" + "=" * 70)
    print("FINDING 4: RY-PATTERN DIFFERENCES = 74")
    print("=" * 70)

    def ry(base):
        return 'Y' if base in 'TC' else 'R'

    patterns = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        pat = ry(codon[0]) + ry(codon[1]) + ry(codon[2])
        patterns[pat].append((codon, aa))

    sums = {}
    for pat in sorted(patterns):
        total = sum(sc(aa) for _, aa in patterns[pat])
        sums[pat] = total
        print(f"  {pat}: sc = {total}  {props(total)}")

    # Complementary pattern pairs
    print("\n  Complementary pattern differences:")
    pairs = [("YYY", "RRR"), ("YYR", "RRY"), ("YRY", "RYR"), ("YRR", "RYY")]
    for p1, p2 in pairs:
        d = sums[p1] - sums[p2]
        s = sums[p1] + sums[p2]
        print(f"  {p1}-{p2} = {sums[p1]}-{sums[p2]} = {d}  {props(abs(d))}")
        print(f"  {p1}+{p2} = {s}  {props(s)}")

    print(f"\n  ★ YYR - RRY = {sums['YYR'] - sums['RRY']} = 74 = weight of constant part!")
    print(f"  ★ This connects the RY-pattern structure to the amino acid constant part.")

    # Are there other pattern differences that equal key numbers?
    print("\n  All pairwise differences:")
    pats = sorted(sums.keys())
    for i, p1 in enumerate(pats):
        for p2 in pats[i + 1:]:
            d = abs(sums[p1] - sums[p2])
            if d % 37 == 0 or is_sm(d):
                print(f"    |{p1}-{p2}| = {d}  {props(d)}")


def finding_5_complement_pairs():
    """FINDING: Reverse complement total = 3404 = 37×92 = 74×46."""
    print("\n" + "=" * 70)
    print("FINDING 5: REVERSE COMPLEMENT STRUCTURE")
    print("=" * 70)

    comp = {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}

    def rev_comp(codon):
        return ''.join(comp[b] for b in reversed(codon))

    # For each codon, pair with its reverse complement
    # Sum sc(codon_aa) + sc(revcomp_aa) for all 32 unique pairs
    seen = set()
    pair_sums = []
    for codon in sorted(GENETIC_CODE):
        rc = rev_comp(codon)
        key = tuple(sorted([codon, rc]))
        if key in seen: continue
        seen.add(key)
        s = sc(GENETIC_CODE[codon]) + sc(GENETIC_CODE[rc])
        pair_sums.append(s)

    total = sum(pair_sums)
    print(f"  Sum of all 32 rev-complement pair sums: {total}")
    print(f"  {props(total)}")
    print(f"  ★ 3404 = 37 × 92 = 74 × 46")
    print(f"  ★ This equals the total of all 64 codons' side chains: "
          f"{sum(sc(aa) for aa in GENETIC_CODE.values())}")

    # It must equal the sum of all 64 codons because each codon appears exactly once
    # in the pair enumeration. So this is just: Σ(sc(all 64 codons)) = 3404 = 37×92
    all_64_sc = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE)
    all_64_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in GENETIC_CODE)
    print(f"\n  Total side chains (all 64): {all_64_sc} = {props(all_64_sc)}")
    print(f"  Total full weight (stop=74, all 64): {all_64_fw74} = {props(all_64_fw74)}")
    print(f"  ★ 8140 = 3700 + 4440 (Octet I + Octet II with stop=74)")
    print(f"  ★ 8140 / 74 = {8140 // 74} = {props(8140 // 74)}")
    print(f"  ★ 8140 / 37 = {8140 // 37} = {props(8140 // 37)}")

    # But also: 3404 = all side chains cell-by-cell
    # 3404 + 64×74 (constant parts, all = 74 including 3 stops at 74) = 3404 + 4736 = 8140
    # Actually: 3404 + N×74 where N = number of non-stop codons = 61
    # plus 3×0 (stop constant parts) = 3404 + 61×74 = 3404 + 4514 = 7918 ≠ 8140
    # Hmm, so the relationship depends on stop assignment.

    # Deeper: can we split the 32 pairs into meaningful groups?
    # Self-complementary codons
    self_comp = []
    for codon in sorted(GENETIC_CODE):
        if codon == rev_comp(codon):
            self_comp.append(codon)
    print(f"\n  Self-complementary codons: {self_comp}")
    if self_comp:
        sc_self = sum(sc(GENETIC_CODE[c]) for c in self_comp)
        print(f"    Their side chain sum: {sc_self} {props(sc_self)}")

    # Pair sum distribution
    print(f"\n  Per-pair sums:")
    seen = set()
    for codon in sorted(GENETIC_CODE):
        rc = rev_comp(codon)
        key = tuple(sorted([codon, rc]))
        if key in seen: continue
        seen.add(key)
        aa1, aa2 = GENETIC_CODE[codon], GENETIC_CODE[rc]
        s = sc(aa1) + sc(aa2)
        p = props(s)
        if p:
            print(f"    {codon}({aa1})↔{rc}({aa2}): {sc(aa1)}+{sc(aa2)}={s} {p}")


def finding_6_all_64_divisibility():
    """The total of all 64 codons' side chains = 3404 = 37×92."""
    print("\n" + "=" * 70)
    print("FINDING 6: TOTAL OF ALL 64 CODONS")
    print("=" * 70)

    total_sc = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE)
    total_fw0 = sum(fw(GENETIC_CODE[c]) for c in GENETIC_CODE)
    total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in GENETIC_CODE)

    print(f"  All 64 codons side chains: {total_sc} = {props(total_sc)}")
    print(f"  All 64 codons full weight (stop=0): {total_fw0} = {props(total_fw0)}")
    print(f"  All 64 codons full weight (stop=74): {total_fw74} = {props(total_fw74)}")
    print(f"  ★ 3404 = 37 × 92 = 4 × 851 = 4 × 23 × 37")
    print(f"  ★ 7918 = 37 × 214 = 2 × 37 × 107")
    print(f"  ★ 8140 = 37 × 220 = 74 × 110 = 20 × 407 = 20 × 11 × 37")
    print(f"  ★ 8140 / 20 = 407 = 11 × 37")
    print(f"  ★ 8140 = 3700 + 4440 (known)")

    # Can we decompose 3404 further?
    # 3404 by first base
    for b in BASES:
        codons = [c for c in GENETIC_CODE if c[0] == b]
        s = sum(sc(GENETIC_CODE[c]) for c in codons)
        print(f"  First={b}: sc={s} {props(s)}")


def finding_7_balanced_splits():
    """Search for exact equalities in natural groupings."""
    print("\n" + "=" * 70)
    print("FINDING 7: EXACT EQUALITIES AND BALANCES")
    print("=" * 70)

    # Already found: 3rd base T = 3rd base C = 864
    # What other equalities exist?

    # All 2-way splits of codons by position + base
    equalities = []

    # Position × specific base
    for pos in range(3):
        sums = {}
        for b in BASES:
            codons = [c for c in GENETIC_CODE if c[pos] == b]
            sums[b] = sum(sc(GENETIC_CODE[c]) for c in codons)
        for i, b1 in enumerate(BASES):
            for b2 in BASES[i + 1:]:
                if sums[b1] == sums[b2]:
                    equalities.append((f"pos{pos + 1}={b1}", f"pos{pos + 1}={b2}", sums[b1]))

    print("  Exact equalities (cell-by-cell side chain sums):")
    for a, b, v in equalities:
        print(f"    ★ {a} = {b} = {v} {props(v)}")

    # Full weight equalities
    print("\n  Full weight equalities (stop=0):")
    for pos in range(3):
        sums = {}
        for b in BASES:
            codons = [c for c in GENETIC_CODE if c[pos] == b]
            sums[b] = sum(fw(GENETIC_CODE[c]) for c in codons)
        for i, b1 in enumerate(BASES):
            for b2 in BASES[i + 1:]:
                if sums[b1] == sums[b2]:
                    print(f"    ★ pos{pos + 1}={b1} = pos{pos + 1}={b2} = {sums[b1]} {props(sums[b1])}")

    # Full weight (stop=74)
    print("\n  Full weight equalities (stop=74):")
    for pos in range(3):
        sums = {}
        for b in BASES:
            codons = [c for c in GENETIC_CODE if c[pos] == b]
            sums[b] = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        for i, b1 in enumerate(BASES):
            for b2 in BASES[i + 1:]:
                if sums[b1] == sums[b2]:
                    print(f"    ★ pos{pos + 1}={b1} = pos{pos + 1}={b2} = {sums[b1]} {props(sums[b1])}")

    # Two-base equalities
    print("\n  Two-position equalities (side chains):")
    for p1 in range(3):
        for p2 in range(p1 + 1, 3):
            # Group by pair of bases at positions p1, p2
            sums = {}
            for b1 in BASES:
                for b2 in BASES:
                    codons = [c for c in GENETIC_CODE if c[p1] == b1 and c[p2] == b2]
                    sums[(b1, b2)] = sum(sc(GENETIC_CODE[c]) for c in codons)
            keys = list(sums.keys())
            for i, k1 in enumerate(keys):
                for k2 in keys[i + 1:]:
                    if sums[k1] == sums[k2] and sums[k1] > 0:
                        print(f"    ★ pos{p1 + 1}{p2 + 1}={k1[0]}{k1[1]} = "
                              f"pos{p1 + 1}{p2 + 1}={k2[0]}{k2[1]} = {sums[k1]} {props(sums[k1])}")


def finding_8_sm_counting_third_base():
    """SM counting applied to third-base groupings."""
    print("\n" + "=" * 70)
    print("FINDING 8: SM COUNTING × THIRD BASE COMBINATIONS")
    print("=" * 70)

    # SM counting by column, but split by third base type
    # For each column, only count the product that appears when third base ∈ Y vs R
    print("  SM counting split by third base Y vs R:")
    for label, b3_set in [("Y (T,C)", "TC"), ("R (A,G)", "AG")]:
        total = 0
        for b1 in BASES:
            for b2 in BASES:
                seen = set()
                for b3 in b3_set:
                    aa = GENETIC_CODE[b1 + b2 + b3]
                    if aa not in seen:
                        seen.add(aa)
                        total += sc(aa)
        print(f"    Third base ∈ {label}: {total} {props(total)}")

    print("\n  SM counting split by third base M vs K:")
    for label, b3_set in [("M (T,G)", "TG"), ("K (C,A)", "CA")]:
        total = 0
        for b1 in BASES:
            for b2 in BASES:
                seen = set()
                for b3 in b3_set:
                    aa = GENETIC_CODE[b1 + b2 + b3]
                    if aa not in seen:
                        seen.add(aa)
                        total += sc(aa)
        print(f"    Third base ∈ {label}: {total} {props(total)}")


def finding_9_triple_37():
    """Look for triple patterns involving 37."""
    print("\n" + "=" * 70)
    print("FINDING 9: THREE-WAY SPLITS GIVING 37-MULTIPLES")
    print("=" * 70)

    # Split 64 codons by first base into 4 groups of 16
    print("  Cell-by-cell side chains by first base:")
    first_sums = {}
    for b in BASES:
        codons = [c for c in GENETIC_CODE if c[0] == b]
        s = sum(sc(GENETIC_CODE[c]) for c in codons)
        first_sums[b] = s
        print(f"    First={b}: {s} {props(s)}")

    # Check all triples (pick 3 of 4 bases)
    print("\n  Three-base combinations (first position):")
    for excluded in BASES:
        included = [b for b in BASES if b != excluded]
        s = sum(first_sums[b] for b in included)
        print(f"    Exclude {excluded}: {s} {props(s)}")

    # Same for second position
    print("\n  Cell-by-cell side chains by second base:")
    second_sums = {}
    for b in BASES:
        codons = [c for c in GENETIC_CODE if c[1] == b]
        s = sum(sc(GENETIC_CODE[c]) for c in codons)
        second_sums[b] = s
        print(f"    Second={b}: {s} {props(s)}")

    # Third position (already done above but repeating for cross-reference)
    print("\n  Cell-by-cell side chains by third base:")
    third_sums = {}
    for b in BASES:
        codons = [c for c in GENETIC_CODE if c[2] == b]
        s = sum(sc(GENETIC_CODE[c]) for c in codons)
        third_sums[b] = s
        print(f"    Third={b}: {s} {props(s)}")

    # Cross-position products
    print("\n  Position sum cross-comparisons:")
    for pos_name, sums in [("1st", first_sums), ("2nd", second_sums), ("3rd", third_sums)]:
        for l1, b1_set in [("Y", "TC"), ("R", "AG"), ("M", "TG"), ("K", "CA")]:
            s = sum(sums[b] for b in b1_set)
            p = props(s)
            if p:
                print(f"    {pos_name} ∈ {l1}: {s} {p}")


def finding_10_hidden_pythagorean():
    """Search for Pythagorean-like relationships."""
    print("\n" + "=" * 70)
    print("FINDING 10: PYTHAGOREAN AND ALGEBRAIC RELATIONSHIPS")
    print("=" * 70)

    # Known: Octet I: 333=37×9, 592=37×16, 925=37×25, and 9+16=25 (not Pythagorean on weights)
    # But 9=3², 16=4², 25=5² and 3²+4²=5² (Egyptian triangle)

    # Can we find similar for other groupings?
    # For any grouping: sc = A, cp = B, fw = A+B. If all three are 37×n, check n values.

    print("  Known: Octet I: sc=333(37×9), cp=592(37×16), fw=925(37×25)")
    print("    9=3², 16=4², 25=5², 3²+4²=5² ✓")

    # Try Octet II (with stop=74 for constant part)
    oi_cols = [(b1, b2) for b1 in BASES for b2 in BASES if is_homogeneous(b1, b2)]
    oii_cols = [(b1, b2) for b1 in BASES for b2 in BASES if not is_homogeneous(b1, b2)]

    oii_prods = sm_products(oii_cols)
    oii_sc = sum(sc(aa) for aa in oii_prods)
    oii_cp0 = sum(0 if aa == 'stop' else 74 for aa in oii_prods)
    oii_cp74 = sum(74 for aa in oii_prods)
    oii_fw0 = oii_sc + oii_cp0
    oii_fw74 = oii_sc + oii_cp74

    print(f"\n  Octet II (stop=0): sc={oii_sc}, cp={oii_cp0}, fw={oii_fw0}")
    print(f"    sc/37={oii_sc / 37:.2f}, cp/37={oii_cp0 / 37:.2f}, fw/37={oii_fw0 / 37:.2f}")
    print(f"  Octet II (stop=74): sc={oii_sc}, cp={oii_cp74}, fw={oii_fw74}")
    print(f"    sc/37={oii_sc / 37:.2f}, cp/37={oii_cp74 / 37:.2f}, fw/37={oii_fw74 / 37:.2f}")

    # Search across multiple groupings
    print("\n  Searching for 37-triple patterns across groupings:")
    groupings = []

    # Per first base
    for b in BASES:
        cols = [(b, b2) for b2 in BASES]
        prods = sm_products(cols)
        sc_val = sum(sc(aa) for aa in prods)
        cp_val = sum(74 for aa in prods if aa != 'stop')
        fw_val = sc_val + cp_val
        groupings.append((f"1st={b}", sc_val, cp_val, fw_val))

    # Per second base
    for b in BASES:
        cols = [(b1, b) for b1 in BASES]
        prods = sm_products(cols)
        sc_val = sum(sc(aa) for aa in prods)
        cp_val = sum(74 for aa in prods if aa != 'stop')
        fw_val = sc_val + cp_val
        groupings.append((f"2nd={b}", sc_val, cp_val, fw_val))

    # Y/R, M/K for each position
    for pos_name, pos in [("1st", 0), ("2nd", 1)]:
        for label, bases in [("Y", "TC"), ("R", "AG"), ("M", "TG"), ("K", "CA")]:
            if pos == 0:
                cols = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in bases]
            else:
                cols = [(b1, b2) for b1 in BASES for b2 in BASES if b2 in bases]
            prods = sm_products(cols)
            sc_val = sum(sc(aa) for aa in prods)
            cp_val = sum(74 for aa in prods if aa != 'stop')
            fw_val = sc_val + cp_val
            groupings.append((f"{pos_name}∈{label}", sc_val, cp_val, fw_val))

    for name, s, c, f in groupings:
        s37 = s % 37 == 0
        c37 = c % 37 == 0
        f37 = f % 37 == 0
        if s37 and c37 and f37:
            a, b, cc = s // 37, c // 37, f // 37
            pyth = "★ PYTHAGOREAN!" if any(
                x * x + y * y == z * z for x, y, z in [(a, b, cc)]
            ) else ""
            print(f"    {name}: sc={s}(37×{a}), cp={c}(37×{b}), fw={f}(37×{cc}) {pyth}")
        elif s37 or f37:
            print(f"    {name}: sc={s}{'(÷37)' if s37 else ''}, cp={c}{'(÷37)' if c37 else ''}, "
                  f"fw={f}{'(÷37)' if f37 else ''}")


def finding_11_exhaustive_sm_search():
    """Exhaustive search for SM numbers in all natural 8-column splits."""
    print("\n" + "=" * 70)
    print("FINDING 11: EXHAUSTIVE SM SEARCH (all 8-column splits)")
    print("=" * 70)

    all_cols = [(b1, b2) for b1 in BASES for b2 in BASES]

    # There are C(16,8) = 12870 ways to split 16 into two groups of 8
    sm_hits = []
    for combo in combinations(range(16), 8):
        g1 = [all_cols[i] for i in combo]
        g2 = [all_cols[i] for i in range(16) if i not in combo]
        p1 = sm_products(g1)
        p2 = sm_products(g2)
        s1 = sum(sc(aa) for aa in p1)
        s2 = sum(sc(aa) for aa in p2)
        both_sm = is_sm(s1) and is_sm(s2)
        both_37 = s1 % 37 == 0 and s2 % 37 == 0
        equal = s1 == s2
        if both_sm or (both_37 and not both_sm) or equal:
            cols1 = [''.join(c) for c in g1]
            cols2 = [''.join(c) for c in g2]
            tag = "BOTH SM" if both_sm else ("BOTH ×37" if both_37 else "EQUAL")
            sm_hits.append((tag, s1, s2, cols1, cols2))

    # Deduplicate (g1/g2 are mirror pairs)
    seen = set()
    for tag, s1, s2, c1, c2 in sm_hits:
        key = (min(s1, s2), max(s1, s2), tag)
        if key in seen: continue
        seen.add(key)
        print(f"\n  {tag}: {s1} ({sm_repr(s1) if is_sm(s1) else props(s1)}) | "
              f"{s2} ({sm_repr(s2) if is_sm(s2) else props(s2)})")
        print(f"    G1: {c1}")
        print(f"    G2: {c2}")

    print(f"\n  Total unique hits: {len(seen)}")


def main():
    print("=" * 70)
    print("DEEP DIVE: HUNTING FOR LEVEL 6")
    print("=" * 70)

    finding_1_third_base_symmetry()
    finding_2_third_base_rumer_sm()
    finding_3_gc_content()
    finding_4_ry_pattern_74()
    finding_5_complement_pairs()
    finding_6_all_64_divisibility()
    finding_7_balanced_splits()
    finding_8_sm_counting_third_base()
    finding_9_triple_37()
    finding_10_hidden_pythagorean()
    finding_11_exhaustive_sm_search()


if __name__ == '__main__':
    main()
