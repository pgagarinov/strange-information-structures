"""
Level 6 Findings: New information structures of the genetic code
not reported by Panov & Filatov (2024) or Shcherbak & Makukov (2013).

These findings were discovered by systematic computational exploration.
"""

from math import gcd
from collections import defaultdict

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

BASES = ['T', 'C', 'A', 'G']


def sc(aa):
    if aa == 'stop': return 0
    W = {'G': 1, 'A': 15, 'S': 31, 'P': 42, 'V': 43, 'T': 45, 'C': 47,
         'L': 57, 'I': 57, 'N': 58, 'D': 59, 'Q': 72, 'K': 72,
         'E': 73, 'M': 75, 'H': 81, 'F': 91, 'R': 100, 'Y': 107, 'W': 130}
    return 41 if aa == 'P' else W[aa]


def fw(aa, stop=0):
    if aa == 'stop': return stop
    FW = {'G': 75, 'A': 89, 'S': 105, 'P': 115, 'V': 117, 'T': 119, 'C': 121,
          'L': 131, 'I': 131, 'N': 132, 'D': 133, 'Q': 146, 'K': 146,
          'E': 147, 'M': 149, 'H': 155, 'F': 165, 'R': 174, 'Y': 181, 'W': 204}
    return FW[aa]


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


def sm_repr(n):
    for m in range(n // 999 + 1):
        rem = n - m * 999
        if rem >= 0 and rem % 111 == 0 and 0 <= rem // 111 <= 9:
            k = rem // 111
            if m == 0: return f"{k} × 111"
            return f"{m} × 999 + {k} × 111" if k > 0 else f"{m} × 999"
    return None


def divider(title):
    print(f"\n{'━' * 70}")
    print(f"  {title}")
    print(f"{'━' * 70}")


def main():
    print("╔" + "═" * 68 + "╗")
    print("║" + " LEVEL 6: NEW INFORMATION STRUCTURES OF THE GENETIC CODE ".center(68) + "║")
    print("║" + " Computational discoveries beyond Panov & Filatov (2024) ".center(68) + "║")
    print("╚" + "═" * 68 + "╝")

    # =========================================================================
    divider("6.1  THE GC-CONTENT SIGNATURE (double SM-number)")
    # =========================================================================
    print("""
  NEW RESOURCE: GC content (number of G or C bases in a codon).
  This is fundamental to molecular biology but never used as a grouping
  principle for weight analysis of the genetic code.

  The 64 codons split by GC content into groups of 8, 24, 24, 8.""")

    gc_groups = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        gc = sum(1 for b in codon if b in 'GC')
        gc_groups[gc].append((codon, aa))

    for gc in sorted(gc_groups):
        codons = gc_groups[gc]
        s = sum(sc(aa) for _, aa in codons)
        f74 = sum(fw(aa, stop=74) for _, aa in codons)
        sm_s = sm_repr(s)
        sm_f = sm_repr(f74)
        marker = " ◄◄◄" if sm_s or sm_f else ""
        print(f"\n    GC = {gc}  ({len(codons):2d} codons):  "
              f"side chains = {s:5d}   full weight(stop=74) = {f74:5d}{marker}")
        if sm_s: print(f"      side chains = {sm_s}  =  37 × {s // 37}" if s % 37 == 0 else f"      side chains = {sm_s}")
        if sm_f: print(f"      full weight = {sm_f}  =  37 × {f74 // 37}" if f74 % 37 == 0 else f"      full weight = {sm_f}")

    print(f"""
  ★ The 24 codons with GC content = 1 produce TWO INDEPENDENT SM-numbers:

        side chains       = 1332 = 999 + 333    (= 37 × 36 = 111 × 12)
        full weight(s=74) = 3108 = 3 × 999 + 111 (= 37 × 84 = 111 × 28)

  ★ The probability of two random integers BOTH being SM-numbers is
    approximately (1/111)² ≈ 8 × 10⁻⁵.

  ★ GC content is one of the most fundamental properties in molecular biology
    (affecting melting temperature, stability, gene expression). That it
    produces SM-numbers suggests the "signal" permeates the code's structure
    more deeply than previously recognized.""")

    # Additional: GC=1 constant parts
    gc1_prods = [aa for _, aa in gc_groups[1]]
    n_stop = sum(1 for aa in gc1_prods if aa == 'stop')
    n_nonstop = len(gc1_prods) - n_stop
    cp = n_nonstop * 74
    print(f"\n    GC=1 has {n_stop} stop codons, {n_nonstop} non-stop codons")
    print(f"    Constant parts (stop=0) = {n_nonstop} × 74 = {cp} = 37 × {cp // 37} = 74 × {cp // 74}")
    print(f"    Constant parts (stop=74) = {len(gc1_prods)} × 74 = {len(gc1_prods) * 74} = {sm_repr(len(gc1_prods) * 74)}")
    print(f"      = 37 × {len(gc1_prods) * 74 // 37}")

    # =========================================================================
    divider("6.2  THIRD-BASE RUMER PAIRS → SM-NUMBERS AND ×74 DIVISIBILITY")
    # =========================================================================
    print("""
  The paper analyzes Rumer pairs M=(T,G) and K=(C,A) for the FIRST and
  SECOND codon bases. The THIRD (wobble) base was never examined.

  Splitting all 64 codons by third base into Rumer pairs:""")

    for label, b3_set in [("M = (T, G)", "TG"), ("K = (C, A)", "CA")]:
        codons = [c for c in GENETIC_CODE if c[2] in b3_set]
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        total_fw0 = sum(fw(GENETIC_CODE[c]) for c in codons)
        total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        print(f"\n    Third base ∈ {label}  (32 codons):")
        print(f"      side chains        = {total_sc:5d}  = {sm_repr(total_sc) or ''}"
              f"  = 37 × {total_sc // 37}" if total_sc % 37 == 0 else
              f"      side chains        = {total_sc}")
        print(f"      full weight (s=0)  = {total_fw0:5d}  = 37 × {total_fw0 // 37}  = 74 × {total_fw0 // 74}")
        print(f"      full weight (s=74) = {total_fw74:5d}  = {sm_repr(total_fw74) or ''}"
              f"  = 37 × {total_fw74 // 37}" if total_fw74 % 37 == 0 else
              f"      full weight (s=74) = {total_fw74}")

    print(f"""
  ★ Third base M: side chains = 1776 = 999 + 777 (SM-number!)
  ★ Third base K: full weight(stop=74) = 3996 = 3 × 999 + 999 (SM-number!)
  ★ ALL SIX values are divisible by both 37 and 74.
  ★ GCD of M and K side chains = {gcd(1776, 1628)} = 4 × 37 = 2 × 74
  ★ The Rumer split yields SM-numbers at ALL THREE codon positions.""")

    # SM counting version
    print("\n  SM-counting version (third base restricted):")
    for label, b3_set in [("M = (T, G)", "TG"), ("K = (C, A)", "CA")]:
        total = 0
        for b1 in BASES:
            for b2 in BASES:
                seen = set()
                for b3 in b3_set:
                    aa = GENETIC_CODE[b1 + b2 + b3]
                    if aa not in seen:
                        seen.add(aa)
                        total += sc(aa)
        sr = sm_repr(total)
        print(f"    Third base ∈ {label}: {total}" +
              (f" = {sr} = 37 × {total // 37}" if sr and total % 37 == 0 else
               f" = {sr}" if sr else ""))

    print(f"""
  ★ SM-counting + third base M = 1443 = 999 + 444 (SM-number!)
     = 37 × 39 = 111 × 13""")

    # =========================================================================
    divider("6.3  THE CONSTANT-PART ECHO: YYR − RRY = 74")
    # =========================================================================

    def ry(base):
        return 'Y' if base in 'TC' else 'R'

    patterns = {}
    for codon, aa in GENETIC_CODE.items():
        pat = ry(codon[0]) + ry(codon[1]) + ry(codon[2])
        patterns.setdefault(pat, []).append(sc(aa))

    sums = {pat: sum(ws) for pat, ws in patterns.items()}

    print("""
  Group all 64 codons by their purine/pyrimidine (R/Y) pattern.
  There are 2³ = 8 patterns, each containing 8 codons.
""")
    for pat in sorted(sums):
        print(f"    {pat}: {sums[pat]:4d}")

    d = sums['YYR'] - sums['RRY']
    print(f"""
  ★ YYR − RRY = {sums['YYR']} − {sums['RRY']} = {d}

  ★ 74 is the weight of the amino acid CONSTANT PART — the number that
    defines the entire framework of the paper's analysis.

  ★ This difference connects the purine/pyrimidine PATTERN structure of
    codons directly to the amino acid constant part, bridging Level 1
    (structural symmetry) and Level 3 (weight analysis).""")

    # =========================================================================
    divider("6.4  FIRST-PYRIMIDINE BALANCE: side chains = constant parts")
    # =========================================================================

    pyr_cols = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in 'TC']
    prods = sm_products(pyr_cols)
    sc_sum = sum(sc(aa) for aa in prods)
    n_nonstop = sum(1 for aa in prods if aa != 'stop')
    cp_sum = n_nonstop * 74

    print(f"""
  SM-counting for all columns with first base = pyrimidine (T or C):

    Side chain sum    = {sc_sum}
    Constant part sum = {n_nonstop} × 74 = {cp_sum}  (stop counted as 0)

  ★ {sc_sum} = {cp_sum} : PERFECT BALANCE""")
    print(f"  ★ = 74 × 11 = 37 × 22")
    print(f"""
  ★ The paper noted this type of balance for Octet II (sc = cp = 1110).
    This is a NEW instance for a DIFFERENT grouping (first-base pyrimidine).

  Verification — the 13 SM-counted products:""")
    for aa in prods:
        s = sc(aa)
        c = 74 if aa != 'stop' else 0
        print(f"    {aa:>4s}:  sc = {s:3d}   cp = {c:2d}")
    print(f"    {'Σ':>4s}:  sc = {sc_sum:3d}   cp = {cp_sum}")

    # =========================================================================
    divider("6.5  THE TRIPLE SM PARTITION")
    # =========================================================================

    all_cols = [(b1, b2) for b1 in BASES for b2 in BASES]
    total_sm = sum(sc(aa) for aa in sm_products(all_cols))

    print(f"""
  The total SM-counted side chain weight across all 16 columns:

    Total = {total_sm} = 999 + 444 = 111 × 13 (SM-number!)

  This total can be partitioned into 8+8 column groups where BOTH halves
  are SM-numbers. Exhaustive search (all C(16,8) = 12,870 splits) finds
  exactly THREE such partitions:

    333 + 1110  =  3×111 + 10×111  (Octet I | Octet II — KNOWN)
    555 +  888  =  5×111 +  8×111  (NEW)
    666 +  777  =  6×111 +  7×111  (NEW)

  ★ These are the ONLY three ways to split the 16 columns into two groups
    of 8 such that both groups give SM-numbers.

  ★ The three partitions tile the SM-number line at 111-spacing:
    3|10, 5|8, 6|7 — note the consecutive pair 6|7!""")

    # Show the columns for 555|888 and 666|777
    from itertools import combinations

    def is_sm(n):
        if n <= 0: return False
        for m in range(n // 999 + 1):
            rem = n - m * 999
            if rem >= 0 and rem % 111 == 0 and rem // 111 <= 9:
                return True
        return False

    counts = defaultdict(int)
    examples = {}
    for combo in combinations(range(16), 8):
        g1 = [all_cols[i] for i in combo]
        g2 = [all_cols[i] for i in range(16) if i not in combo]
        s1 = sum(sc(aa) for aa in sm_products(g1))
        s2 = sum(sc(aa) for aa in sm_products(g2))
        if is_sm(s1) and is_sm(s2) and s1 <= s2:
            key = (s1, s2)
            counts[key] += 1
            if key not in examples:
                examples[key] = ([''.join(c) for c in g1], [''.join(c) for c in g2])

    for (s1, s2), n in sorted(counts.items()):
        c1, c2 = examples[(s1, s2)]
        known = " (= Octet I | Octet II)" if s1 == 333 else ""
        print(f"\n    {sm_repr(s1):>17s} | {sm_repr(s2):<17s}  "
              f"({n} column grouping{'s' if n > 1 else ''}){known}")
        print(f"    Example: {c1}")
        print(f"             {c2}")

    # =========================================================================
    # 888|555: explain the Rumer-pair swap
    # =========================================================================
    print(f"""
  The 888|555 partition has a beautiful explanation:

    Start from the first-base Y/R split (from the paper):
      Y-first side chains = 814 = 37 × 22
      R-first side chains = 629 = 37 × 17

    Now SWAP one Rumer pair of columns: CT ↔ AG
      CT products (SM): L, sc = 57
      AG products (SM): S + R, sc = 131

    After the swap:
      814 − 57 + 131 = 888 = 8 × 111 (SM!)
      629 + 57 − 131 = 555 = 5 × 111 (SM!)

  ★ The swap cost is 131 − 57 = 74 = the constant part weight!
  ★ One Rumer pair swap, costing exactly 74, converts a ×37 split into a ×111 split.""")

    # =========================================================================
    divider("6.6  UNIVERSAL ×74 DIVISIBILITY OF THE FULL CODE")
    # =========================================================================

    total_sc = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE)
    total_fw0 = sum(fw(GENETIC_CODE[c]) for c in GENETIC_CODE)
    total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in GENETIC_CODE)

    print(f"""
  Summing over all 64 codons (cell-by-cell counting):

    Side chains       = {total_sc} = 37 × {total_sc // 37} = 74 × {total_sc // 74}
    Full weight(s=0)  = {total_fw0} = 37 × {total_fw0 // 37} = 74 × {total_fw0 // 74}
    Full weight(s=74) = {total_fw74} = 37 × {total_fw74 // 37} = 74 × {total_fw74 // 74}

  ★ 8140 / 20 = 407 = 11 × 37 (where 20 = number of amino acids)
  ★ 8140 = 3700 + 4440 (Octet I + Octet II with stop=74, known)
  ★ 3404 = 37 × 92 = 74 × 46

  The number 37 permeates EVERY level of the code's weight structure.""")

    # =========================================================================
    divider("SUMMARY")
    # =========================================================================

    print("""
  Six new information signatures, organized by the resources they require:

  6.1  GC-content grouping → DOUBLE SM-number (1332, 3108)
       New resource: base composition (GC content)
       Independent of codon position

  6.2  Third-base Rumer pairs → SM-numbers + universal ×74
       Extension of paper's M/K analysis to the wobble position
       1776 = 999+777, 3996 = 3×999+999

  6.3  YYR − RRY = 74 (constant part weight)
       Bridges structural symmetry and weight analysis

  6.4  First-pyrimidine sc = cp = 814 = 74 × 11
       New instance of the side-chain/constant-part balance

  6.5  Triple SM partition: 333|1110, 555|888, 666|777
       Only 3 ways to split 16 columns into equal SM-number groups
       The 888|555 split costs exactly ONE Rumer pair = 74

  6.6  All 64 codons: side chains = 3404 = 74 × 46 = 37 × 92
       The number 37 governs the entire code, not just special subsets

  These findings suggest that the information structures extend beyond the
  five levels identified by Panov & Filatov, encompassing base composition
  (GC content), the wobble position, and global properties of the code.""")


if __name__ == '__main__':
    main()
