"""
Hunt for Level 6: Undiscovered information structures in the genetic code.

Systematically explore groupings of codons/amino acids that Panov & Filatov
did not examine, looking for SM-numbers, divisibility by 37, unexpected
equalities, and other patterns.
"""

from math import gcd
from itertools import combinations, product
from collections import defaultdict

# === DATA (from verify.py) ===

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
    """Side chain with proline activation key."""
    if aa == 'stop': return 0
    return 41 if aa == 'P' else AMINO_ACIDS[aa]


def fw(aa, stop=0):
    if aa == 'stop': return stop
    return FULL_WEIGHT[aa]


def is_homogeneous(b1, b2):
    products = [GENETIC_CODE[b1 + b2 + b3] for b3 in BASES]
    return len(set(products)) == 1


def is_sm_number(n):
    """Check if n = m*999 + k*111 for non-negative m, k with 0 <= k <= 9."""
    if n < 0: return False
    for m in range(n // 999 + 1):
        rem = n - m * 999
        if rem >= 0 and rem % 111 == 0:
            k = rem // 111
            if 0 <= k <= 9:
                return True
    return False


def sm_repr(n):
    """Return SM representation string, or None."""
    if n < 0: return None
    for m in range(n // 999 + 1):
        rem = n - m * 999
        if rem >= 0 and rem % 111 == 0:
            k = rem // 111
            if 0 <= k <= 9:
                if m == 0:
                    return f"{k}×111"
                return f"{m}×999+{k}×111" if k > 0 else f"{m}×999"
    return None


def interesting(n):
    """Check if a number has interesting properties."""
    props = []
    if n == 0:
        return []
    sr = sm_repr(n)
    if sr:
        props.append(f"SM({sr})")
    if n > 0 and n % 37 == 0:
        props.append(f"37×{n // 37}")
    if n > 0 and n % 74 == 0:
        props.append(f"74×{n // 74}")
    if n > 0 and n % 111 == 0:
        props.append(f"111×{n // 111}")
    # Check repdigit
    s = str(abs(n))
    if len(s) >= 3 and len(set(s)) == 1:
        props.append(f"repdigit")
    # Digit sum divisible by 9
    if n > 0:
        ds = sum(int(d) for d in str(n))
        if ds % 9 == 0 and n % 9 == 0:
            props.append(f"9×{n // 9}")
    return props


def sm_products(columns):
    """SM counting: each distinct product per column counted once."""
    products = []
    for b1, b2 in columns:
        seen = set()
        for b3 in BASES:
            aa = GENETIC_CODE[b1 + b2 + b3]
            if aa not in seen:
                seen.add(aa)
                products.append(aa)
    return products


# =============================================================================
# EXPLORATION 1: Third-base groupings (unexplored by the paper)
# =============================================================================

def explore_third_base():
    print("\n" + "=" * 70)
    print("EXPLORATION 1: Third-base groupings")
    print("=" * 70)
    print("The paper groups by 1st/2nd base. What about the 3rd (wobble) base?")

    # For each third base, sum side chains of all codons with that third base
    for b3 in BASES:
        codons = [c for c in GENETIC_CODE if c[2] == b3]
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        total_fw = sum(fw(GENETIC_CODE[c]) for c in codons)
        total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        props_sc = interesting(total_sc)
        props_fw = interesting(total_fw)
        props_fw74 = interesting(total_fw74)
        print(f"\n  Third base = {b3} ({len(codons)} codons):")
        print(f"    Side chains (cell-by-cell): {total_sc} {props_sc}")
        print(f"    Full weight (stop=0):       {total_fw} {props_fw}")
        print(f"    Full weight (stop=74):      {total_fw74} {props_fw74}")

    # Third base Y vs R
    for label, bases in [("Y (T,C)", "TC"), ("R (A,G)", "AG")]:
        codons = [c for c in GENETIC_CODE if c[2] in bases]
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        print(f"\n  Third base ∈ {label} ({len(codons)} codons):")
        print(f"    Side chains: {total_sc} {interesting(total_sc)}")
        print(f"    Full weight (stop=74): {total_fw74} {interesting(total_fw74)}")

    # Third base M vs K (Rumer pairs)
    for label, bases in [("M (T,G)", "TG"), ("K (C,A)", "CA")]:
        codons = [c for c in GENETIC_CODE if c[2] in bases]
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        print(f"\n  Third base ∈ {label} ({len(codons)} codons):")
        print(f"    Side chains: {total_sc} {interesting(total_sc)}")
        print(f"    Full weight (stop=74): {total_fw74} {interesting(total_fw74)}")


# =============================================================================
# EXPLORATION 2: Full RY-pattern of triplet
# =============================================================================

def explore_ry_patterns():
    print("\n" + "=" * 70)
    print("EXPLORATION 2: Codon RY-pattern (purine/pyrimidine triplet pattern)")
    print("=" * 70)

    def ry(base):
        return 'Y' if base in 'TC' else 'R'

    patterns = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        pat = ry(codon[0]) + ry(codon[1]) + ry(codon[2])
        patterns[pat].append((codon, aa))

    for pat in sorted(patterns):
        codons = patterns[pat]
        total_sc = sum(sc(aa) for _, aa in codons)
        total_fw74 = sum(fw(aa, stop=74) for _, aa in codons)
        aas = sorted(set(aa for _, aa in codons if aa != 'stop'))
        props = interesting(total_sc)
        print(f"\n  Pattern {pat} ({len(codons)} codons): {', '.join(aas)}")
        print(f"    Side chains (cell-by-cell): {total_sc} {props}")
        print(f"    Full weight (stop=74):      {total_fw74} {interesting(total_fw74)}")

    # Check YYY+RRR vs YRY+RYR etc (complementary patterns)
    print("\n  --- Complementary pattern pairs ---")
    for p1, p2 in [("YYY", "RRR"), ("YYR", "RRY"), ("YRY", "RYR"), ("YRR", "RYY")]:
        s1 = sum(sc(aa) for _, aa in patterns[p1])
        s2 = sum(sc(aa) for _, aa in patterns[p2])
        props = interesting(s1 + s2)
        eq = interesting(abs(s1 - s2))
        print(f"  {p1}={s1} + {p2}={s2} = {s1 + s2} {props}; diff={s1 - s2} {eq}")


# =============================================================================
# EXPLORATION 3: GC content
# =============================================================================

def explore_gc_content():
    print("\n" + "=" * 70)
    print("EXPLORATION 3: GC content of codons")
    print("=" * 70)

    gc_groups = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        gc = sum(1 for b in codon if b in 'GC')
        gc_groups[gc].append((codon, aa))

    for gc in sorted(gc_groups):
        codons = gc_groups[gc]
        total_sc = sum(sc(aa) for _, aa in codons)
        total_fw74 = sum(fw(aa, stop=74) for _, aa in codons)
        print(f"\n  GC content = {gc} ({len(codons)} codons):")
        print(f"    Side chains: {total_sc} {interesting(total_sc)}")
        print(f"    Full weight: {total_fw74} {interesting(total_fw74)}")

    # AT-rich (gc=0,1) vs GC-rich (gc=2,3)
    at_rich = sum(sc(aa) for gc in [0, 1] for _, aa in gc_groups[gc])
    gc_rich = sum(sc(aa) for gc in [2, 3] for _, aa in gc_groups[gc])
    print(f"\n  AT-rich (gc≤1): {at_rich} {interesting(at_rich)}")
    print(f"  GC-rich (gc≥2): {gc_rich} {interesting(gc_rich)}")
    print(f"  Sum: {at_rich + gc_rich} {interesting(at_rich + gc_rich)}")
    print(f"  Diff: {at_rich - gc_rich} {interesting(abs(at_rich - gc_rich))}")


# =============================================================================
# EXPLORATION 4: Degeneracy classes
# =============================================================================

def explore_degeneracy():
    print("\n" + "=" * 70)
    print("EXPLORATION 4: Amino acids grouped by degeneracy (codon count)")
    print("=" * 70)

    codon_count = defaultdict(int)
    for aa in GENETIC_CODE.values():
        if aa != 'stop':
            codon_count[aa] += 1

    deg_groups = defaultdict(list)
    for aa, count in codon_count.items():
        deg_groups[count].append(aa)

    for deg in sorted(deg_groups):
        aas = sorted(deg_groups[deg], key=lambda a: AMINO_ACIDS[a])
        total_sc = sum(sc(aa) for aa in aas)
        total_fw = sum(fw(aa) for aa in aas)
        n = len(aas)
        print(f"\n  {deg}-fold degenerate ({n} amino acids): {', '.join(aas)}")
        print(f"    Side chains: {total_sc} {interesting(total_sc)}")
        print(f"    Full weight: {total_fw} {interesting(total_fw)}")
        if n > 0:
            print(f"    Mean side chain: {total_sc / n:.1f}")

    # Special: 1-fold (M, W) vs 6-fold (L, R, S)
    one_fold = sum(sc(aa) for aa in deg_groups[1])
    six_fold = sum(sc(aa) for aa in deg_groups[6])
    print(f"\n  1-fold sum: {one_fold}, 6-fold sum: {six_fold}")
    print(f"  1+6 fold: {one_fold + six_fold} {interesting(one_fold + six_fold)}")


# =============================================================================
# EXPLORATION 5: Complementary codon pairs
# =============================================================================

def explore_complement():
    print("\n" + "=" * 70)
    print("EXPLORATION 5: Reverse complement codon pairs")
    print("=" * 70)

    comp = {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}

    def rev_comp(codon):
        return ''.join(comp[b] for b in reversed(codon))

    pairs_seen = set()
    pair_diffs = []
    pair_sums = []
    for codon in sorted(GENETIC_CODE):
        rc = rev_comp(codon)
        pair_key = tuple(sorted([codon, rc]))
        if pair_key in pairs_seen:
            continue
        pairs_seen.add(pair_key)
        aa1 = GENETIC_CODE[codon]
        aa2 = GENETIC_CODE[rc]
        s1, s2 = sc(aa1), sc(aa2)
        pair_sums.append(s1 + s2)
        pair_diffs.append(abs(s1 - s2))

    total_sum = sum(pair_sums)
    print(f"  Total sum of all pair sums: {total_sum} {interesting(total_sum)}")
    print(f"  Total sum of all |diffs|:   {sum(pair_diffs)} {interesting(sum(pair_diffs))}")

    # Self-complementary codons (codon = rev_comp)
    print("\n  Self-complementary codons:")
    for codon in sorted(GENETIC_CODE):
        if codon == rev_comp(codon):
            aa = GENETIC_CODE[codon]
            print(f"    {codon} → {aa} (sc={sc(aa)})")

    # Sum amino acid weights encoded by complementary codon pairs
    # Group: for each column, pair it with its Rumer-complement column
    print("\n  Complementary column pairs (SM counting):")
    done = set()
    for b1 in BASES:
        for b2 in BASES:
            col = (b1, b2)
            rc_col = (comp[b2], comp[b1])  # reverse complement of column
            pair_key = tuple(sorted([col, rc_col]))
            if pair_key in done:
                continue
            done.add(pair_key)
            prods1 = set(GENETIC_CODE[b1 + b2 + b3] for b3 in BASES)
            prods2 = set(GENETIC_CODE[rc_col[0] + rc_col[1] + b3] for b3 in BASES)
            s1 = sum(sc(aa) for aa in prods1)
            s2 = sum(sc(aa) for aa in prods2)
            if col == rc_col:
                print(f"    {''.join(col)} (self): {s1} {interesting(s1)}")
            else:
                print(f"    {''.join(col)}↔{''.join(rc_col)}: {s1}+{s2}={s1 + s2} {interesting(s1 + s2)}, diff={s1 - s2}")


# =============================================================================
# EXPLORATION 6: Diagonal patterns in the 4×4 table
# =============================================================================

def explore_diagonals():
    print("\n" + "=" * 70)
    print("EXPLORATION 6: Diagonal and geometric patterns in the code table")
    print("=" * 70)

    # Main diagonal: b1 == b2
    diag_main = [(b, b) for b in BASES]
    prods = sm_products(diag_main)
    total = sum(sc(aa) for aa in prods)
    print(f"\n  Main diagonal (b1=b2): columns {[''.join(c) for c in diag_main]}")
    print(f"    Side chains (SM): {total} {interesting(total)}")

    # Anti-diagonal: b1+b2 indices sum to 3
    anti_diag = [(BASES[i], BASES[3 - i]) for i in range(4)]
    prods = sm_products(anti_diag)
    total = sum(sc(aa) for aa in prods)
    print(f"\n  Anti-diagonal: columns {[''.join(c) for c in anti_diag]}")
    print(f"    Side chains (SM): {total} {interesting(total)}")

    # Complement diagonal: (b1, comp(b1))
    comp = {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}
    comp_diag = [(b, comp[b]) for b in BASES]
    prods = sm_products(comp_diag)
    total = sum(sc(aa) for aa in prods)
    print(f"\n  Complement diagonal (b2=comp(b1)): columns {[''.join(c) for c in comp_diag]}")
    print(f"    Side chains (SM): {total} {interesting(total)}")

    # Rumer diagonal: (b1, rumer(b1))
    rumer = {'T': 'G', 'G': 'T', 'C': 'A', 'A': 'C'}
    rumer_diag = [(b, rumer[b]) for b in BASES]
    prods = sm_products(rumer_diag)
    total = sum(sc(aa) for aa in prods)
    print(f"\n  Rumer diagonal (b2=rumer(b1)): columns {[''.join(c) for c in rumer_diag]}")
    print(f"    Side chains (SM): {total} {interesting(total)}")


# =============================================================================
# EXPLORATION 7: Individual base sums (SM counting, per base position)
# =============================================================================

def explore_per_base():
    print("\n" + "=" * 70)
    print("EXPLORATION 7: Per-base side chain sums (SM counting)")
    print("=" * 70)

    for pos_name, pos in [("first", 0), ("second", 1)]:
        print(f"\n  --- Grouped by {pos_name} base ---")
        for b in BASES:
            if pos == 0:
                cols = [(b, b2) for b2 in BASES]
            else:
                cols = [(b1, b) for b1 in BASES]
            prods = sm_products(cols)
            total = sum(sc(aa) for aa in prods)
            print(f"    {pos_name}={b}: {total} {interesting(total)}")

    # All four 2-base type combinations for Y/R
    print("\n  --- 2x2 grid: first Y/R × second Y/R ---")
    for f_label, f_bases in [("Y", "TC"), ("R", "AG")]:
        for s_label, s_bases in [("Y", "TC"), ("R", "AG")]:
            cols = [(b1, b2) for b1 in BASES for b2 in BASES
                    if b1 in f_bases and b2 in s_bases]
            prods = sm_products(cols)
            total = sum(sc(aa) for aa in prods)
            print(f"    1st∈{f_label}, 2nd∈{s_label}: {total} {interesting(total)}")

    # Cross: M/K for one position, Y/R for the other
    print("\n  --- Cross: first M/K × second Y/R ---")
    for f_label, f_bases in [("M", "TG"), ("K", "CA")]:
        for s_label, s_bases in [("Y", "TC"), ("R", "AG")]:
            cols = [(b1, b2) for b1 in BASES for b2 in BASES
                    if b1 in f_bases and b2 in s_bases]
            prods = sm_products(cols)
            total = sum(sc(aa) for aa in prods)
            print(f"    1st∈{f_label}, 2nd∈{s_label}: {total} {interesting(total)}")

    print("\n  --- Cross: first Y/R × second M/K ---")
    for f_label, f_bases in [("Y", "TC"), ("R", "AG")]:
        for s_label, s_bases in [("M", "TG"), ("K", "CA")]:
            cols = [(b1, b2) for b1 in BASES for b2 in BASES
                    if b1 in f_bases and b2 in s_bases]
            prods = sm_products(cols)
            total = sum(sc(aa) for aa in prods)
            print(f"    1st∈{f_label}, 2nd∈{s_label}: {total} {interesting(total)}")


# =============================================================================
# EXPLORATION 8: Modular arithmetic patterns
# =============================================================================

def explore_modular():
    print("\n" + "=" * 70)
    print("EXPLORATION 8: Modular arithmetic of amino acid side chains")
    print("=" * 70)

    weights = sorted(set(sc(aa) for aa in AMINO_ACIDS))
    print(f"\n  Unique side chain weights (activated): {weights}")
    print(f"  Sum of all 20: {sum(sc(aa) for aa in AMINO_ACIDS)} "
          f"{interesting(sum(sc(aa) for aa in AMINO_ACIDS))}")

    # Residues mod 37
    print(f"\n  Weights mod 37:")
    for aa in sorted(AMINO_ACIDS, key=lambda a: AMINO_ACIDS[a]):
        w = sc(aa)
        print(f"    {aa}: {w} ≡ {w % 37} (mod 37)")

    residues = [sc(aa) % 37 for aa in AMINO_ACIDS]
    print(f"\n  Sum of residues mod 37: {sum(residues)} ≡ {sum(residues) % 37} (mod 37)")
    print(f"  Product of residues mod 37: {1}")  # too large, skip

    # Digital roots
    print(f"\n  Digital roots of side chain weights:")
    dr_groups = defaultdict(list)
    for aa in sorted(AMINO_ACIDS, key=lambda a: AMINO_ACIDS[a]):
        w = sc(aa)
        dr = w % 9 if w % 9 != 0 else (9 if w > 0 else 0)
        dr_groups[dr].append(aa)
    for dr in sorted(dr_groups):
        aas = dr_groups[dr]
        total = sum(sc(aa) for aa in aas)
        print(f"    DR={dr}: {', '.join(aas)} → sum={total} {interesting(total)}")


# =============================================================================
# EXPLORATION 9: Nucleotide base weights interaction
# =============================================================================

def explore_nucleotide_weights():
    print("\n" + "=" * 70)
    print("EXPLORATION 9: Nucleotide base molecular weights")
    print("=" * 70)
    print("  Base weights: T=126, C=111, A=135, G=151")

    BASE_WT = {'T': 126, 'C': 111, 'A': 135, 'G': 151}

    # C=111 is itself an SM number!
    print(f"\n  C weight = 111 = 1×111 (SM number!)")
    print(f"  T weight = 126 = {126 // 9}×9")
    print(f"  A weight = 135 = {135 // 9}×9")
    print(f"  G weight = 151 (prime)")
    print(f"  T+C+A+G = {sum(BASE_WT.values())} {interesting(sum(BASE_WT.values()))}")
    print(f"  Y sum (T+C) = {BASE_WT['T'] + BASE_WT['C']}")
    print(f"  R sum (A+G) = {BASE_WT['A'] + BASE_WT['G']}")
    print(f"  M sum (T+G) = {BASE_WT['T'] + BASE_WT['G']}")
    print(f"  K sum (C+A) = {BASE_WT['C'] + BASE_WT['A']}")

    # Codon nucleotide weight vs amino acid weight
    print("\n  Codon nucleotide weight vs amino acid side chain:")
    for b1 in BASES:
        for b2 in BASES:
            cols = [(b1, b2)]
            prods = sm_products(cols)
            nuc_wt = BASE_WT[b1] + BASE_WT[b2]
            aa_sc = sum(sc(aa) for aa in prods)
            ratio = aa_sc / nuc_wt if nuc_wt else 0
            print(f"    {b1}{b2}: nuc_wt={nuc_wt}, sc_sum={aa_sc}, ratio={ratio:.3f}")


# =============================================================================
# EXPLORATION 10: Exhaustive search for SM-number partitions
# =============================================================================

def explore_partition_search():
    print("\n" + "=" * 70)
    print("EXPLORATION 10: Systematic search for equal-sum partitions")
    print("=" * 70)

    all_cols = [(b1, b2) for b1 in BASES for b2 in BASES]

    # For each way to split the 16 columns into 2 groups of 8,
    # check if both sums are SM numbers (too many to enumerate all C(16,8)=12870)
    # Instead, check natural splits

    splits = {
        "Octet I/II (known)": (
            [(b1, b2) for b1, b2 in all_cols if is_homogeneous(b1, b2)],
            [(b1, b2) for b1, b2 in all_cols if not is_homogeneous(b1, b2)]
        ),
        "1st base Y/R": (
            [(b1, b2) for b1, b2 in all_cols if b1 in 'TC'],
            [(b1, b2) for b1, b2 in all_cols if b1 in 'AG']
        ),
        "2nd base Y/R": (
            [(b1, b2) for b1, b2 in all_cols if b2 in 'TC'],
            [(b1, b2) for b1, b2 in all_cols if b2 in 'AG']
        ),
        "1st base M/K (known)": (
            [(b1, b2) for b1, b2 in all_cols if b1 in 'TG'],
            [(b1, b2) for b1, b2 in all_cols if b1 in 'CA']
        ),
        "2nd base M/K (known)": (
            [(b1, b2) for b1, b2 in all_cols if b2 in 'TG'],
            [(b1, b2) for b1, b2 in all_cols if b2 in 'CA']
        ),
        "Complement pair (b1,comp(b1)) vs rest": (
            [(b, {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}[b]) for b in BASES]
            + [(b, b) for b in BASES],
            [(b1, b2) for b1, b2 in all_cols
             if (b1, b2) not in set([(b, {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}[b]) for b in BASES]
                                    + [(b, b) for b in BASES])]
        ),
        "Main+anti diag vs off-diag": (
            [(BASES[i], BASES[i]) for i in range(4)] + [(BASES[i], BASES[3 - i]) for i in range(4)],
            [(b1, b2) for b1, b2 in all_cols
             if (b1, b2) not in set([(BASES[i], BASES[i]) for i in range(4)]
                                    + [(BASES[i], BASES[3 - i]) for i in range(4)])]
        ),
    }

    for name, (g1, g2) in splits.items():
        p1 = sm_products(g1)
        p2 = sm_products(g2)
        s1 = sum(sc(aa) for aa in p1)
        s2 = sum(sc(aa) for aa in p2)
        i1 = interesting(s1)
        i2 = interesting(s2)
        si = interesting(s1 + s2)
        eq = "EQUAL!" if s1 == s2 else ""
        marker = " ★" if (i1 or i2 or s1 == s2) else ""
        print(f"\n  {name}:{marker}")
        print(f"    Group 1 ({len(g1)} cols): {s1} {i1}")
        print(f"    Group 2 ({len(g2)} cols): {s2} {i2}")
        print(f"    Sum: {s1 + s2} {si}  Diff: {abs(s1 - s2)} {interesting(abs(s1 - s2))}")
        if eq:
            print(f"    {eq}")


# =============================================================================
# EXPLORATION 11: Three-base combinations beyond pairs
# =============================================================================

def explore_three_bases():
    print("\n" + "=" * 70)
    print("EXPLORATION 11: Codons with exactly 3 distinct bases")
    print("=" * 70)

    # 24 codons with all 3 bases different
    all_diff = [c for c in GENETIC_CODE if len(set(c)) == 3]
    all_same = [c for c in GENETIC_CODE if len(set(c)) == 1]
    two_same = [c for c in GENETIC_CODE if len(set(c)) == 2
                and sum(1 for a, b in [(c[0], c[1]), (c[0], c[2]), (c[1], c[2])] if a == b) == 1]

    for label, codons in [("All 3 different", all_diff),
                          ("All 3 same", all_same),
                          ("Exactly 2 same (known)", two_same)]:
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        total_fw74 = sum(fw(GENETIC_CODE[c], stop=74) for c in codons)
        print(f"\n  {label} ({len(codons)} codons):")
        print(f"    Side chains: {total_sc} {interesting(total_sc)}")
        print(f"    Full weight (stop=74): {total_fw74} {interesting(total_fw74)}")

    # Among all-different: group by which base is missing
    print("\n  All-3-different, grouped by missing base:")
    for missing in BASES:
        present = [b for b in BASES if b != missing]
        codons = [c for c in all_diff if missing not in c]
        # Actually: codons that use exactly the 3 bases in 'present'
        codons = [c for c in all_diff if set(c) == set(present)]
        total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
        print(f"    Missing {missing}: {len(codons)} codons, sc={total_sc} {interesting(total_sc)}")

    # Among all-different: group by which base is the singleton (appears once)
    # In a 3-different codon, each base appears exactly once, so this = all permutations
    # Group by whether the codon has more purines or pyrimidines
    print("\n  All-3-different, by purine count:")
    for n_pur in range(4):
        codons = [c for c in all_diff if sum(1 for b in c if b in 'AG') == n_pur]
        if codons:
            total_sc = sum(sc(GENETIC_CODE[c]) for c in codons)
            print(f"    {n_pur} purines: {len(codons)} codons, sc={total_sc} {interesting(total_sc)}")


# =============================================================================
# EXPLORATION 12: Amino acid chemical groupings
# =============================================================================

def explore_chemical_groups():
    print("\n" + "=" * 70)
    print("EXPLORATION 12: Chemical property groupings")
    print("=" * 70)

    # Standard biochemical groupings
    groups = {
        "Nonpolar aliphatic": ['G', 'A', 'V', 'L', 'I', 'P', 'M'],
        "Aromatic": ['F', 'W', 'Y'],
        "Polar uncharged": ['S', 'T', 'C', 'N', 'Q'],
        "Positively charged": ['K', 'R', 'H'],
        "Negatively charged": ['D', 'E'],
        # Alternative groupings
        "Sulfur-containing": ['C', 'M'],
        "Hydroxyl-containing": ['S', 'T', 'Y'],
        "Amide-containing": ['N', 'Q'],
        "Branched-chain": ['V', 'L', 'I'],
        "Small (sc≤15)": ['G', 'A'],
        "Tiny (sc≤31)": ['G', 'A', 'S'],
    }

    for name, aas in groups.items():
        total_sc = sum(sc(aa) for aa in aas)
        total_fw = sum(fw(aa) for aa in aas)
        print(f"\n  {name}: {', '.join(aas)}")
        print(f"    Side chains: {total_sc} {interesting(total_sc)}")
        print(f"    Full weight: {total_fw} {interesting(total_fw)}")


# =============================================================================
# EXPLORATION 13: Codon value encoding → weight correlations
# =============================================================================

def explore_codon_encoding():
    print("\n" + "=" * 70)
    print("EXPLORATION 13: Base position × weight interactions")
    print("=" * 70)

    # Assign numerical values to bases
    encodings = {
        "T=0,C=1,A=2,G=3": {'T': 0, 'C': 1, 'A': 2, 'G': 3},
        "T=1,C=2,A=3,G=4": {'T': 1, 'C': 2, 'A': 3, 'G': 4},
        "Rumer: T=0,G=0,C=1,A=1": {'T': 0, 'G': 0, 'C': 1, 'A': 1},
        "Y/R: T=0,C=0,A=1,G=1": {'T': 0, 'C': 0, 'A': 1, 'G': 1},
    }

    for enc_name, enc in encodings.items():
        print(f"\n  Encoding: {enc_name}")
        total = 0
        for codon, aa in GENETIC_CODE.items():
            codon_val = sum(enc[b] * (4 ** (2 - i)) for i, b in enumerate(codon))
            weight = sc(aa)
            total += codon_val * weight
        print(f"    Σ(codon_val × side_chain) = {total} {interesting(total)}")

        # Weighted by position
        for pos in range(3):
            total_pos = 0
            for codon, aa in GENETIC_CODE.items():
                total_pos += enc[codon[pos]] * sc(aa)
            print(f"    Σ(base[{pos}]_val × side_chain) = {total_pos} {interesting(total_pos)}")


# =============================================================================
# EXPLORATION 14: Amino acid pairs and symmetries
# =============================================================================

def explore_aa_pairs():
    print("\n" + "=" * 70)
    print("EXPLORATION 14: Amino acid pairing symmetries")
    print("=" * 70)

    # Pairs that share a column (differ only in 3rd base wobble)
    print("\n  Amino acid pairs sharing a column:")
    col_pairs = set()
    for b1 in BASES:
        for b2 in BASES:
            aas_in_col = sorted(set(GENETIC_CODE[b1 + b2 + b3] for b3 in BASES) - {'stop'})
            if len(aas_in_col) == 2:
                pair = tuple(aas_in_col)
                s = sc(pair[0]) + sc(pair[1])
                d = abs(sc(pair[0]) - sc(pair[1]))
                col_pairs.add(pair)
                print(f"    {b1}{b2}: {pair[0]}({sc(pair[0])}) + {pair[1]}({sc(pair[1])}) "
                      f"= {s} {interesting(s)}, diff={d} {interesting(d)}")

    # Weight of each amino acid paired with its Rumer-transformed partner
    rumer = {'T': 'G', 'G': 'T', 'C': 'A', 'A': 'C'}
    print("\n  Rumer partner amino acid pairs (same column position):")
    seen = set()
    for b1 in BASES:
        for b2 in BASES:
            rb1, rb2 = rumer[b1], rumer[b2]
            key = tuple(sorted([(b1, b2), (rb1, rb2)]))
            if key in seen: continue
            seen.add(key)
            prods1 = sm_products([(b1, b2)])
            prods2 = sm_products([(rb1, rb2)])
            s1 = sum(sc(aa) for aa in prods1)
            s2 = sum(sc(aa) for aa in prods2)
            print(f"    {b1}{b2}↔{rb1}{rb2}: {s1}+{s2}={s1 + s2} {interesting(s1 + s2)}")


# =============================================================================
# EXPLORATION 15: Weight palindromes and digit patterns
# =============================================================================

def explore_digit_patterns():
    print("\n" + "=" * 70)
    print("EXPLORATION 15: Weight digit patterns and concatenation")
    print("=" * 70)

    # Concatenate all side chain weights in genetic code table order
    print("\n  Side chain weights in table order (TCAG × TCAG):")
    row_sums = {}
    for b1 in BASES:
        weights = []
        for b2 in BASES:
            prods = sm_products([(b1, b2)])
            ws = [sc(aa) for aa in prods]
            weights.extend(ws)
        row_sum = sum(weights)
        row_sums[b1] = row_sum
        print(f"    {b1}-row: sum={row_sum} {interesting(row_sum)}")

    print(f"\n  Row sums: {dict(row_sums)}")
    for b1, b2 in [('T', 'G'), ('C', 'A')]:
        s = row_sums[b1] + row_sums[b2]
        print(f"    {b1}+{b2} = {s} {interesting(s)}")

    # Sum of squares of side chains
    sq_sum = sum(sc(aa) ** 2 for aa in AMINO_ACIDS)
    print(f"\n  Sum of squares of side chains: {sq_sum} {interesting(sq_sum)}")

    # Product patterns
    prod_all = 1
    for aa in AMINO_ACIDS:
        prod_all *= sc(aa)
    print(f"  Product of all side chains: {prod_all}")
    print(f"    = {prod_all} {interesting(prod_all)}")


# =============================================================================
# EXPLORATION 16: Full weight SM counting — all natural 2-way splits
# =============================================================================

def explore_full_weight_sm():
    print("\n" + "=" * 70)
    print("EXPLORATION 16: Full weight sums (SM counting, stop=74)")
    print("=" * 70)

    all_cols = [(b1, b2) for b1 in BASES for b2 in BASES]

    splits = [
        ("1st base Y vs R", 'TC', 'AG', 0),
        ("1st base M vs K", 'TG', 'CA', 0),
        ("2nd base Y vs R", 'TC', 'AG', 1),
        ("2nd base M vs K", 'TG', 'CA', 1),
    ]

    for name, g1_bases, g2_bases, pos in splits:
        if pos == 0:
            cols1 = [(b1, b2) for b1, b2 in all_cols if b1 in g1_bases]
            cols2 = [(b1, b2) for b1, b2 in all_cols if b1 in g2_bases]
        else:
            cols1 = [(b1, b2) for b1, b2 in all_cols if b2 in g1_bases]
            cols2 = [(b1, b2) for b1, b2 in all_cols if b2 in g2_bases]
        s1 = sum(fw(aa, stop=74) for aa in sm_products(cols1))
        s2 = sum(fw(aa, stop=74) for aa in sm_products(cols2))
        print(f"\n  {name}:")
        print(f"    {g1_bases}: fw={s1} {interesting(s1)}")
        print(f"    {g2_bases}: fw={s2} {interesting(s2)}")
        print(f"    Sum: {s1 + s2} {interesting(s1 + s2)}  Diff: {abs(s1 - s2)} {interesting(abs(s1 - s2))}")
        print(f"    GCD: {gcd(s1, s2)} {interesting(gcd(s1, s2))}")


# =============================================================================
# EXPLORATION 17: "Virtual" amino acid 21 = stop
# =============================================================================

def explore_stop_variations():
    print("\n" + "=" * 70)
    print("EXPLORATION 17: What stop weight makes more SM-numbers?")
    print("=" * 70)

    oii_cols = [(b1, b2) for b1 in BASES for b2 in BASES if not is_homogeneous(b1, b2)]

    for stop_sc in range(0, 150):
        def test_sc(aa):
            if aa == 'stop': return stop_sc
            return 41 if aa == 'P' else AMINO_ACIDS[aa]

        oii_prods = sm_products(oii_cols)
        s = sum(test_sc(aa) for aa in oii_prods)
        if is_sm_number(s) and stop_sc != 0:
            print(f"  stop_sc={stop_sc}: Octet II side chains = {s} ({sm_repr(s)})")


# =============================================================================
# EXPLORATION 18: Amino acid number assignments
# =============================================================================

def explore_numbering():
    print("\n" + "=" * 70)
    print("EXPLORATION 18: Natural numbering of amino acids")
    print("=" * 70)

    # Number amino acids 1-20 by weight
    sorted_aas = sorted(AMINO_ACIDS, key=lambda a: AMINO_ACIDS[a])
    print("  Amino acids by side chain weight:")
    for i, aa in enumerate(sorted_aas, 1):
        print(f"    {i:2d}. {aa} (sc={sc(aa)})")

    # Sum of numbers for octet I vs octet II amino acids
    oi_cols = [(b1, b2) for b1 in BASES for b2 in BASES if is_homogeneous(b1, b2)]
    oii_cols = [(b1, b2) for b1 in BASES for b2 in BASES if not is_homogeneous(b1, b2)]

    numbering = {aa: i + 1 for i, aa in enumerate(sorted_aas)}

    oi_aas = set(aa for aa in sm_products(oi_cols) if aa != 'stop')
    oii_aas = set(aa for aa in sm_products(oii_cols) if aa != 'stop')

    oi_nums = sum(numbering[aa] for aa in oi_aas)
    oii_nums = sum(numbering[aa] for aa in oii_aas)
    print(f"\n  Octet I amino acids: {sorted(oi_aas, key=lambda a: AMINO_ACIDS[a])}")
    print(f"    Sum of numbers: {oi_nums} {interesting(oi_nums)}")
    print(f"  Octet II amino acids: {sorted(oii_aas, key=lambda a: AMINO_ACIDS[a])}")
    print(f"    Sum of numbers: {oii_nums} {interesting(oii_nums)}")
    print(f"    Total: {oi_nums + oii_nums} {interesting(oi_nums + oii_nums)}")

    # Number squared
    oi_sq = sum(numbering[aa] ** 2 for aa in oi_aas)
    oii_sq = sum(numbering[aa] ** 2 for aa in oii_aas)
    print(f"\n  Sum of squares: Octet I={oi_sq}, Octet II={oii_sq}")
    print(f"    Octet I: {interesting(oi_sq)}")
    print(f"    Octet II: {interesting(oii_sq)}")


# =============================================================================
# EXPLORATION 19: Three-way and four-way splits
# =============================================================================

def explore_multiway_splits():
    print("\n" + "=" * 70)
    print("EXPLORATION 19: Three-way and four-way column splits")
    print("=" * 70)

    all_cols = [(b1, b2) for b1 in BASES for b2 in BASES]

    # Four-way by first base
    print("\n  Four-way by first base (SM side chains):")
    for b in BASES:
        cols = [(b, b2) for b2 in BASES]
        total = sum(sc(aa) for aa in sm_products(cols))
        print(f"    First={b}: {total} {interesting(total)}")

    # Four-way by second base
    print("\n  Four-way by second base:")
    for b in BASES:
        cols = [(b1, b) for b1 in BASES]
        total = sum(sc(aa) for aa in sm_products(cols))
        print(f"    Second={b}: {total} {interesting(total)}")

    # Three-way: choose 3 of 4 bases for first position
    print("\n  Three-way: first base in {3 of 4 bases}:")
    for excluded in BASES:
        bases = [b for b in BASES if b != excluded]
        cols = [(b1, b2) for b1 in bases for b2 in BASES]
        total = sum(sc(aa) for aa in sm_products(cols))
        print(f"    Exclude {excluded}: {total} {interesting(total)}")


# =============================================================================
# EXPLORATION 20: Sum over unique amino acids (no SM counting)
# =============================================================================

def explore_unique_aa():
    print("\n" + "=" * 70)
    print("EXPLORATION 20: Unique amino acid sums (each AA counted once)")
    print("=" * 70)

    all_20_sc = sum(sc(aa) for aa in AMINO_ACIDS)
    all_20_fw = sum(fw(aa) for aa in AMINO_ACIDS)
    print(f"  All 20 AAs side chains: {all_20_sc} {interesting(all_20_sc)}")
    print(f"  All 20 AAs full weight: {all_20_fw} {interesting(all_20_fw)}")
    print(f"  All 20 AAs full weight + 20×74 constant = {all_20_fw}")
    print(f"  Side chains + 20×74 = {all_20_sc + 20 * 74} vs full = {all_20_fw}")

    # Known: all_20_sc with proline key
    # What about without proline key?
    all_20_sc_native = sum(AMINO_ACIDS[aa] for aa in AMINO_ACIDS)
    print(f"\n  All 20 AAs side chains (native, no proline key): {all_20_sc_native} {interesting(all_20_sc_native)}")

    # Octet-specific unique AAs
    oi_aas = set()
    oii_aas = set()
    for b1 in BASES:
        for b2 in BASES:
            for b3 in BASES:
                aa = GENETIC_CODE[b1 + b2 + b3]
                if aa == 'stop': continue
                if is_homogeneous(b1, b2):
                    oi_aas.add(aa)
                else:
                    oii_aas.add(aa)

    shared = oi_aas & oii_aas
    only_i = oi_aas - oii_aas
    only_ii = oii_aas - oi_aas
    print(f"\n  Octet I only: {sorted(only_i)} sc={sum(sc(aa) for aa in only_i)} {interesting(sum(sc(aa) for aa in only_i))}")
    print(f"  Octet II only: {sorted(only_ii)} sc={sum(sc(aa) for aa in only_ii)} {interesting(sum(sc(aa) for aa in only_ii))}")
    print(f"  Shared: {sorted(shared)} sc={sum(sc(aa) for aa in shared)} {interesting(sum(sc(aa) for aa in shared))}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("HUNTING FOR LEVEL 6")
    print("Exploring undiscovered information structures of the genetic code")
    print("=" * 70)

    explore_third_base()
    explore_ry_patterns()
    explore_gc_content()
    explore_degeneracy()
    explore_complement()
    explore_diagonals()
    explore_per_base()
    explore_modular()
    explore_nucleotide_weights()
    explore_partition_search()
    explore_three_bases()
    explore_chemical_groups()
    explore_codon_encoding()
    explore_aa_pairs()
    explore_digit_patterns()
    explore_full_weight_sm()
    explore_stop_variations()
    explore_numbering()
    explore_multiway_splits()
    explore_unique_aa()


if __name__ == '__main__':
    main()
