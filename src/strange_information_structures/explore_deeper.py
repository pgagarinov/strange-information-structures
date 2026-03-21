"""
Deeper explorations for Level 6+.
"""

from math import gcd
from collections import defaultdict
from itertools import combinations
import sys

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
AA_SC = {
    'G': 1, 'A': 15, 'S': 31, 'P': 42, 'V': 43, 'T': 45, 'C': 47,
    'L': 57, 'I': 57, 'N': 58, 'D': 59, 'Q': 72, 'K': 72,
    'E': 73, 'M': 75, 'H': 81, 'F': 91, 'R': 100, 'Y': 107, 'W': 130,
}
AA_FW = {
    'G': 75, 'A': 89, 'S': 105, 'P': 115, 'V': 117, 'T': 119, 'C': 121,
    'L': 131, 'I': 131, 'N': 132, 'D': 133, 'Q': 146, 'K': 146,
    'E': 147, 'M': 149, 'H': 155, 'F': 165, 'R': 174, 'Y': 181, 'W': 204,
}
ALL_AA = sorted(AA_SC.keys(), key=lambda a: AA_SC[a])


def sc(aa):
    if aa == 'stop': return 0
    return 41 if aa == 'P' else AA_SC[aa]


def fw(aa, stop=0):
    return stop if aa == 'stop' else AA_FW[aa]


def is_sm(n):
    if n <= 0: return False
    for m in range(n // 999 + 1):
        rem = n - m * 999
        if rem >= 0 and rem % 111 == 0 and rem // 111 <= 9:
            return True
    return False


def sm_repr(n):
    if n <= 0: return None
    for m in range(n // 999 + 1):
        rem = n - m * 999
        if rem >= 0 and rem % 111 == 0 and 0 <= rem // 111 <= 9:
            k = rem // 111
            if m == 0: return f"{k}×111"
            return f"{m}×999+{k}×111" if k > 0 else f"{m}×999"
    return None


def props(n):
    parts = []
    sr = sm_repr(n)
    if sr: parts.append(f"SM({sr})")
    if n > 0 and n % 111 == 0: parts.append(f"111×{n // 111}")
    elif n > 0 and n % 74 == 0: parts.append(f"74×{n // 74}")
    elif n > 0 and n % 37 == 0: parts.append(f"37×{n // 37}")
    return ' '.join(parts)


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


def divider(title):
    print(f"\n{'━' * 70}")
    print(f"  {title}")
    print(f"{'━' * 70}")


# =============================================================================
# A: The 1776 coincidence
# =============================================================================

def explore_1776():
    divider("A. THE 1776 COINCIDENCE")
    print("""
  Two completely different calculations yield the same number:

  (1) Paper Sig 3.5: First-pyrimidine full weight (SM counting, stop=74) = 1776
  (2) New Finding 6.2: Third-base M=(T,G) side chains (cell-by-cell) = 1776
  (3) New Finding 6.1: GC=1 constant parts (stop=74) = 24×74 = 1776

  All three = 999 + 777 = 37 × 48 = 74 × 24 = 111 × 16
""")
    # Verify all three
    # (1) Paper Sig 3.5
    pyr_cols = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in 'TC']
    pyr_prods = sm_products(pyr_cols)
    v1 = sum(fw(aa, stop=74) for aa in pyr_prods)
    print(f"  (1) First-pyr full weight (SM, stop=74): {v1}")

    # (2) Third base M
    v2 = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE if c[2] in 'TG')
    print(f"  (2) Third-base M side chains (cell):     {v2}")

    # (3) GC=1 constant parts with stop=74
    gc1 = [c for c in GENETIC_CODE if sum(1 for b in c if b in 'GC') == 1]
    v3 = len(gc1) * 74
    print(f"  (3) GC=1 constant parts (stop=74):       {v3}  ({len(gc1)} × 74)")

    print(f"\n  All equal 1776: {v1 == v2 == v3 == 1776}")

    # Is there a structural reason?
    print(f"""
  ★ Three routes to 1776 from three independent principles:
    - Position-specific (first base pyrimidine) × full weight × SM counting
    - Position-specific (third base Rumer-M) × side chain × cell-by-cell
    - Position-free (GC content = 1) × constant part × cell-by-cell

  ★ 1776 = 24 × 74 = 48 × 37 = 16 × 111

  ★ The number 1776 is also the year of American independence.
    (Coincidence, but memorable.)""")


# =============================================================================
# B: 4×4 matrix properties
# =============================================================================

def explore_matrix():
    divider("B. THE 4×4 SIDE-CHAIN MATRIX")

    # Build the matrix M[b1][b2] = SM-counted side chain sum for column (b1,b2)
    M = {}
    for b1 in BASES:
        for b2 in BASES:
            prods = sm_products([(b1, b2)])
            M[(b1, b2)] = sum(sc(aa) for aa in prods)

    print("\n  SM-counted side chain sums per column:\n")
    print("        T     C     A     G   | Row Σ")
    print("       ─── ───── ───── ─────  ├──────")
    row_sums = {}
    for b1 in BASES:
        vals = [M[(b1, b2)] for b2 in BASES]
        row_sums[b1] = sum(vals)
        print(f"  {b1}:  {vals[0]:>4d}  {vals[1]:>4d}  {vals[2]:>4d}  {vals[3]:>4d}  | {sum(vals):>4d}  {props(sum(vals))}")

    col_sums = {}
    print("       ─── ───── ───── ─────  ├──────")
    cs = []
    for b2 in BASES:
        s = sum(M[(b1, b2)] for b1 in BASES)
        col_sums[b2] = s
        cs.append(s)
    print(f"  Σ:  {cs[0]:>4d}  {cs[1]:>4d}  {cs[2]:>4d}  {cs[3]:>4d}  | {sum(cs):>4d}")
    print(f"       {props(cs[0]):>4s}  {props(cs[1]):>4s}  {props(cs[2]):>4s}  {props(cs[3]):>4s}")

    # Diagonal sums
    main_diag = sum(M[(BASES[i], BASES[i])] for i in range(4))
    anti_diag = sum(M[(BASES[i], BASES[3 - i])] for i in range(4))
    print(f"\n  Main diagonal (TT+CC+AA+GG): {main_diag}  {props(main_diag)}")
    print(f"  Anti-diagonal (TG+CA+AC+GT): {anti_diag}  {props(anti_diag)}")

    # Rumer pairs: sum of M[x,y] + M[R(x),R(y)]
    rumer = {'T': 'G', 'G': 'T', 'C': 'A', 'A': 'C'}
    print(f"\n  Rumer-paired column sums M[x,y] + M[R(x),R(y)]:")
    seen = set()
    for b1 in BASES:
        for b2 in BASES:
            rb1, rb2 = rumer[b1], rumer[b2]
            key = tuple(sorted([(b1, b2), (rb1, rb2)]))
            if key in seen: continue
            seen.add(key)
            s = M[(b1, b2)] + M[(rb1, rb2)]
            print(f"    {b1}{b2}+{rb1}{rb2}: {M[(b1, b2)]}+{M[(rb1, rb2)]}={s}  {props(s)}")

    # Symmetry: M[b1,b2] vs M[b2,b1] (transpose)
    print(f"\n  Transpose symmetry M[x,y] vs M[y,x]:")
    for i, b1 in enumerate(BASES):
        for b2 in BASES[i + 1:]:
            d = M[(b1, b2)] - M[(b2, b1)]
            if d != 0:
                print(f"    {b1}{b2}={M[(b1, b2)]}, {b2}{b1}={M[(b2, b1)]}, diff={d}  {props(abs(d))}")

    # Row/column Rumer pair sums
    print(f"\n  Row Rumer pair sums:")
    for b1, rb1 in [('T', 'G'), ('C', 'A')]:
        s = row_sums[b1] + row_sums[rb1]
        print(f"    Row {b1} + Row {rb1}: {row_sums[b1]}+{row_sums[rb1]}={s}  {props(s)}")


# =============================================================================
# C: Amino acid subset search for SM numbers
# =============================================================================

def explore_aa_subsets():
    divider("C. AMINO ACID SUBSETS SUMMING TO SM-NUMBERS")

    weights = {aa: sc(aa) for aa in ALL_AA}
    total = sum(weights.values())
    print(f"\n  Total of all 20 side chains: {total}  {props(total)}")

    # Check all subsets of size k for SM-number sums
    print(f"\n  Searching subsets of size k where sum is SM-number...")
    print(f"  (Showing only subsets with interesting structure)\n")

    sm_subsets = defaultdict(list)
    for k in range(2, 19):
        for combo in combinations(ALL_AA, k):
            s = sum(weights[aa] for aa in combo)
            if is_sm(s):
                sm_subsets[k].append((s, combo))

    for k in sorted(sm_subsets):
        entries = sm_subsets[k]
        if len(entries) > 30:
            # Just show SM-number values and counts
            val_counts = defaultdict(int)
            for s, _ in entries:
                val_counts[s] += 1
            print(f"  k={k}: {len(entries)} subsets → SM values: "
                  f"{', '.join(f'{sm_repr(v)}({n})' for v, n in sorted(val_counts.items()))}")
        else:
            for s, combo in sorted(entries):
                print(f"  k={k}: {list(combo)} = {s} = {sm_repr(s)}")

    # Special: subsets that equal exactly 333, 666, 777, 888, 999
    print(f"\n  Subsets summing to key SM-numbers:")
    for target in [111, 222, 333, 444, 555, 666, 777, 888, 999]:
        found = []
        for k in range(1, 21):
            for combo in combinations(ALL_AA, k):
                if sum(weights[aa] for aa in combo) == target:
                    found.append(combo)
        if found:
            print(f"\n    {target} = {sm_repr(target)}:")
            for combo in found[:5]:  # limit output
                print(f"      {list(combo)} (k={len(combo)})")
            if len(found) > 5:
                print(f"      ... and {len(found) - 5} more")


# =============================================================================
# D: Symmetrized code (TGA: stop → C)
# =============================================================================

def explore_symmetrized():
    divider("D. THE SYMMETRIZED CODE (TGA: stop → Cys)")

    sym_code = dict(GENETIC_CODE)
    sym_code['TGA'] = 'C'

    print("  The paper discusses the symmetrized code (Section 6) where")
    print("  TGA encodes Cys instead of stop. This exists in nature")
    print("  (euplotid nuclear code). Checking weight signatures:\n")

    # Octets don't change (TG column was heterogeneous and stays heterogeneous)
    oi = [(b1, b2) for b1 in BASES for b2 in BASES
          if len(set(sym_code[b1 + b2 + b3] for b3 in BASES)) == 1]
    oii = [(b1, b2) for b1 in BASES for b2 in BASES
           if len(set(sym_code[b1 + b2 + b3] for b3 in BASES)) > 1]

    print(f"  Octet I columns: {len(oi)} ({'same' if len(oi) == 8 else 'CHANGED!'})")
    print(f"  Octet II columns: {len(oii)}")

    def sym_sc(aa):
        if aa == 'stop': return 0
        return 41 if aa == 'P' else AA_SC[aa]

    def sym_sm_products(columns):
        products = []
        for b1, b2 in columns:
            seen = set()
            for b3 in BASES:
                aa = sym_code[b1 + b2 + b3]
                if aa not in seen:
                    seen.add(aa)
                    products.append(aa)
        return products

    # Key signatures
    oi_prods = sym_sm_products(oi)
    oii_prods = sym_sm_products(oii)
    oi_sc = sum(sym_sc(aa) for aa in oi_prods)
    oii_sc = sum(sym_sc(aa) for aa in oii_prods)

    print(f"\n  Octet I side chains (SM): {oi_sc}  {props(oi_sc)}")
    print(f"  Octet II side chains (SM): {oii_sc}  {props(oii_sc)}")

    # Cell-by-cell
    oi_cell = sum(sym_sc(sym_code[b1 + b2 + b3]) for b1, b2 in oi for b3 in BASES)
    oii_cell = sum(sym_sc(sym_code[b1 + b2 + b3]) for b1, b2 in oii for b3 in BASES)
    oi_fw = sum(fw(sym_code[b1 + b2 + b3]) for b1, b2 in oi for b3 in BASES)
    oii_fw_74 = sum(fw(sym_code[b1 + b2 + b3], stop=74) for b1, b2 in oii for b3 in BASES)
    oii_fw_0 = sum(fw(sym_code[b1 + b2 + b3], stop=0) for b1, b2 in oii for b3 in BASES)

    print(f"\n  Cell-by-cell:")
    print(f"    Octet I full weight: {oi_fw}  {props(oi_fw)}")
    print(f"    Octet II full weight (stop=0): {oii_fw_0}  {props(oii_fw_0)}")
    print(f"    Octet II full weight (stop=74): {oii_fw_74}  {props(oii_fw_74)}")
    print(f"    GCD(I, II stop=74): {gcd(oi_fw, oii_fw_74)}  {props(gcd(oi_fw, oii_fw_74))}")

    # GC=1 in symmetrized code
    gc1 = [(c, sym_code[c]) for c in sym_code if sum(1 for b in c if b in 'GC') == 1]
    gc1_sc = sum(sym_sc(aa) for _, aa in gc1)
    gc1_fw = sum(fw(aa, stop=74) for _, aa in gc1)
    print(f"\n  GC=1 in symmetrized code:")
    print(f"    side chains: {gc1_sc}  {props(gc1_sc)}")
    print(f"    full weight (stop=74): {gc1_fw}  {props(gc1_fw)}")

    # Total
    total_sc = sum(sym_sc(sym_code[c]) for c in sym_code)
    total_fw = sum(fw(sym_code[c], stop=74) for c in sym_code)
    print(f"\n  All 64 codons:")
    print(f"    side chains: {total_sc}  {props(total_sc)}")
    print(f"    full weight: {total_fw}  {props(total_fw)}")

    # New signatures unique to symmetrized code?
    # 36 codons with 2 identical bases
    print(f"\n  Signature 3.2 in symmetrized code:")
    pyrimidine_pair = []
    purine_pair = []
    for codon in sym_code:
        b1, b2, b3 = codon
        if sum(a == b for a, b in [(b1, b2), (b1, b3), (b2, b3)]) == 1:
            pair_base = b1 if b1 == b2 or b1 == b3 else b2
            if pair_base in 'TC':
                pyrimidine_pair.append(sym_sc(sym_code[codon]))
            else:
                purine_pair.append(sym_sc(sym_code[codon]))
    print(f"    Pyrimidine pair sum: {sum(pyrimidine_pair)}  {props(sum(pyrimidine_pair))}")
    print(f"    Purine pair sum: {sum(purine_pair)}  {props(sum(purine_pair))}")


# =============================================================================
# E: Amino acid PAIRS summing to key numbers
# =============================================================================

def explore_aa_pairs():
    divider("E. AMINO ACID PAIRS AND COMPLEMENTS")

    print("\n  Pairs of amino acids whose side chains sum to 111:")
    for i, a1 in enumerate(ALL_AA):
        for a2 in ALL_AA[i + 1:]:
            s = sc(a1) + sc(a2)
            if s == 111:
                print(f"    {a1}({sc(a1)}) + {a2}({sc(a2)}) = 111")

    print("\n  Pairs summing to 74:")
    for i, a1 in enumerate(ALL_AA):
        for a2 in ALL_AA[i + 1:]:
            s = sc(a1) + sc(a2)
            if s == 74:
                print(f"    {a1}({sc(a1)}) + {a2}({sc(a2)}) = 74")

    print("\n  Pairs summing to 37:")
    for i, a1 in enumerate(ALL_AA):
        for a2 in ALL_AA[i + 1:]:
            if sc(a1) + sc(a2) == 37:
                print(f"    {a1}({sc(a1)}) + {a2}({sc(a2)}) = 37")

    # Anticodon pairing weights
    comp = {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}

    print("\n  Codon vs anticodon amino acid weight differences:")
    print("  (anticodon = reverse complement, read 3'→5')")
    total_same = 0
    total_diff = 0
    for b1 in BASES:
        for b2 in BASES:
            # Column anticodon: reverse complement of (b1,b2,*) reading
            # Anticodon of b1-b2-b3 is comp(b3)-comp(b2)-comp(b1)
            # The first two bases of the anticodon column are comp(b2), comp(b1)
            anti_col = (comp[b2], comp[b1])
            prods_orig = sm_products([(b1, b2)])
            prods_anti = sm_products([anti_col])
            s_orig = sum(sc(aa) for aa in prods_orig)
            s_anti = sum(sc(aa) for aa in prods_anti)
            d = s_orig - s_anti
            if abs(d) % 37 == 0 and d != 0:
                print(f"    {b1}{b2} (sc={s_orig}) ↔ anti {''.join(anti_col)} (sc={s_anti}): "
                      f"diff={d} {props(abs(d))}")


# =============================================================================
# F: Position-weighted sums
# =============================================================================

def explore_position_weighted():
    divider("F. POSITION-WEIGHTED SUMS")

    print("  Weight amino acid by its codon's base position values.")
    print("  Assign T=0, C=1, A=2, G=3.\n")

    base_val = {'T': 0, 'C': 1, 'A': 2, 'G': 3}

    for pos in range(3):
        total = sum(base_val[c[pos]] * sc(GENETIC_CODE[c]) for c in GENETIC_CODE)
        print(f"  Σ(base[{pos + 1}] × side_chain): {total}  {props(total)}")

    # Weighted by complementary value
    comp_val = {'T': 2, 'C': 3, 'A': 0, 'G': 1}  # complement's value
    for pos in range(3):
        total = sum(comp_val[c[pos]] * sc(GENETIC_CODE[c]) for c in GENETIC_CODE)
        print(f"  Σ(comp_val[{pos + 1}] × side_chain): {total}  {props(total)}")

    # Binary: Y=0, R=1
    print()
    for pos in range(3):
        total_y = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE if c[pos] in 'TC')
        total_r = sum(sc(GENETIC_CODE[c]) for c in GENETIC_CODE if c[pos] in 'AG')
        print(f"  Position {pos + 1}: Y-sum={total_y} {props(total_y)}, "
              f"R-sum={total_r} {props(total_r)}, diff={total_y - total_r} {props(abs(total_y - total_r))}")


# =============================================================================
# G: Three-base interaction sums
# =============================================================================

def explore_three_base_interactions():
    divider("G. THREE-POSITION INTERACTION PATTERNS")

    # For each of the 64 codons, assign a "signature" based on
    # whether each position is M or K (Rumer pair)
    print("  Codons grouped by M/K pattern at all three positions:")
    print("  M={T,G}, K={C,A}\n")

    def mk(base):
        return 'M' if base in 'TG' else 'K'

    mk_groups = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        pat = mk(codon[0]) + mk(codon[1]) + mk(codon[2])
        mk_groups[pat].append(sc(aa))

    for pat in sorted(mk_groups):
        vals = mk_groups[pat]
        s = sum(vals)
        print(f"    {pat} ({len(vals):2d} codons): sc = {s:5d}  {props(s)}")

    # Complementary MK pairs
    opp = {'M': 'K', 'K': 'M'}
    print("\n  Complementary MK-pattern pairs:")
    seen = set()
    for pat in sorted(mk_groups):
        anti = ''.join(opp[c] for c in pat)
        key = tuple(sorted([pat, anti]))
        if key in seen: continue
        seen.add(key)
        s1 = sum(mk_groups[pat])
        s2 = sum(mk_groups[anti])
        d = abs(s1 - s2)
        print(f"    {pat}({s1}) + {anti}({s2}) = {s1 + s2} {props(s1 + s2)}, diff={d} {props(d)}")


# =============================================================================
# H: Deeper GC-content analysis
# =============================================================================

def explore_gc_deeper():
    divider("H. DEEPER GC-CONTENT PATTERNS")

    print("  GC content × Octet membership:\n")
    for gc in range(4):
        for octet_name, test in [("I", lambda b1, b2: len(set(GENETIC_CODE[b1 + b2 + b3] for b3 in BASES)) == 1),
                                  ("II", lambda b1, b2: len(set(GENETIC_CODE[b1 + b2 + b3] for b3 in BASES)) > 1)]:
            codons = []
            for c in GENETIC_CODE:
                if sum(1 for b in c if b in 'GC') == gc:
                    b1, b2 = c[0], c[1]
                    if test(b1, b2):
                        codons.append(c)
            if codons:
                total = sum(sc(GENETIC_CODE[c]) for c in codons)
                print(f"    GC={gc}, Octet {octet_name}: {len(codons):2d} codons, sc={total:5d}  {props(total)}")

    # GC content × Y/R first base
    print("\n  GC content × first base type:\n")
    for gc in range(4):
        for label, bases in [("Y-first", "TC"), ("R-first", "AG")]:
            codons = [c for c in GENETIC_CODE
                      if sum(1 for b in c if b in 'GC') == gc and c[0] in bases]
            if codons:
                total = sum(sc(GENETIC_CODE[c]) for c in codons)
                print(f"    GC={gc}, {label}: {len(codons):2d} codons, sc={total:5d}  {props(total)}")

    # GC content × M/K third base
    print("\n  GC content × third base Rumer pair:\n")
    for gc in range(4):
        for label, bases in [("3rd-M", "TG"), ("3rd-K", "CA")]:
            codons = [c for c in GENETIC_CODE
                      if sum(1 for b in c if b in 'GC') == gc and c[2] in bases]
            if codons:
                total = sum(sc(GENETIC_CODE[c]) for c in codons)
                print(f"    GC={gc}, {label}: {len(codons):2d} codons, sc={total:5d}  {props(total)}")


# =============================================================================
# I: Difference patterns between sorted amino acids
# =============================================================================

def explore_differences():
    divider("I. CONSECUTIVE DIFFERENCES AND PARTIAL SUMS")

    sorted_aa = sorted(ALL_AA, key=lambda a: sc(a))
    weights = [sc(aa) for aa in sorted_aa]

    print("  Side chains sorted:")
    print(f"  {list(zip(sorted_aa, weights))}\n")

    diffs = [weights[i + 1] - weights[i] for i in range(len(weights) - 1)]
    print(f"  Consecutive differences: {diffs}")
    print(f"  Sum of differences: {sum(diffs)} (= max - min = {weights[-1]} - {weights[0]})")

    # Running partial sums
    print(f"\n  Partial sums from smallest:")
    partial = 0
    for i, (aa, w) in enumerate(zip(sorted_aa, weights)):
        partial += w
        p = props(partial)
        marker = " ◄" if p else ""
        print(f"    Σ[1..{i + 1:2d}] = {partial:5d}  (+{aa}={w}){marker}  {p}")


# =============================================================================
# J: Column permutation invariants
# =============================================================================

def explore_permutation_invariants():
    divider("J. PERMUTATION-INVARIANT PROPERTIES")

    print("  For each column, compute side chain sum, then sort columns by sum.")
    print("  This is invariant under permutation of bases on the table sides.\n")

    col_sums = []
    for b1 in BASES:
        for b2 in BASES:
            prods = sm_products([(b1, b2)])
            s = sum(sc(aa) for aa in prods)
            col_sums.append((s, b1 + b2))

    col_sums.sort()
    print("  Sorted column SM-counted side chains:")
    for s, col in col_sums:
        print(f"    {col}: {s}")

    vals = [s for s, _ in col_sums]
    print(f"\n  As a sequence: {vals}")
    print(f"  Sum: {sum(vals)}  {props(sum(vals))}")
    print(f"  Min+Max: {vals[0]}+{vals[-1]} = {vals[0] + vals[-1]}  {props(vals[0] + vals[-1])}")

    # Pairs: smallest+largest, second+second-to-last, etc
    print(f"\n  Mirror pairing (1st+16th, 2nd+15th, ...):")
    for i in range(8):
        s = vals[i] + vals[15 - i]
        p = props(s)
        marker = " ◄" if p else ""
        print(f"    {vals[i]:>4d} + {vals[15 - i]:>4d} = {s:>4d}{marker}  {p}")


def main():
    explore_1776()
    explore_matrix()
    explore_aa_subsets()
    explore_symmetrized()
    explore_aa_pairs()
    explore_position_weighted()
    explore_three_base_interactions()
    explore_gc_deeper()
    explore_differences()
    explore_permutation_invariants()


if __name__ == '__main__':
    main()
