"""
Verify all numerical claims from Panov & Filatov (2024):
"Are the Strange Information Structures of the Genetic Code an Accident or an Artifact?"
"""

from math import gcd
from itertools import permutations

# =============================================================================
# DATA FOUNDATION
# =============================================================================

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

# letter -> (name, side_chain_weight, full_weight)
AMINO_ACIDS = {
    'G': ('Glycine',        1,  75),
    'A': ('Alanine',       15,  89),
    'S': ('Serine',        31, 105),
    'P': ('Proline',       42, 115),
    'V': ('Valine',        43, 117),
    'T': ('Threonine',     45, 119),
    'C': ('Cysteine',      47, 121),
    'L': ('Leucine',       57, 131),
    'I': ('Isoleucine',    57, 131),
    'N': ('Asparagine',    58, 132),
    'D': ('Aspartic acid', 59, 133),
    'Q': ('Glutamine',     72, 146),
    'K': ('Lysine',        72, 146),
    'E': ('Glutamic acid', 73, 147),
    'M': ('Methionine',    75, 149),
    'H': ('Histidine',     81, 155),
    'F': ('Phenylalanine', 91, 165),
    'R': ('Arginine',     100, 174),
    'Y': ('Tyrosine',     107, 181),
    'W': ('Tryptophan',   130, 204),
}

# Chain link β weights (Table 20 for Octet I, Table 21 for Octet II)
BETA = {
    'G': 1, 'A': 15, 'S': 14, 'P': 27, 'V': 13, 'T': 13, 'L': 14, 'R': 14,
    'F': 14, 'Y': 14, 'C': 14, 'W': 14, 'H': 14, 'Q': 14, 'I': 13, 'M': 14,
    'N': 14, 'K': 14, 'D': 14, 'E': 14, 'stop': 0,
}

# γ weights from Table 21 (Octet II context)
GAMMA_OCTET_II = {
    'F': 12, 'L': 13, 'Y': 12, 'C': 33, 'W': 13, 'H': 12,
    'Q': 14, 'I': 29, 'M': 14, 'N': 12, 'K': 14, 'S': 17,
    'R': 14, 'D': 12, 'E': 14, 'stop': 0,
}

# β+γ values from Table 22 (first pyrimidine context) — differs from Table 21
# for W (γ=12 vs 13) and adds P (γ=14)
BETA_GAMMA_FIRST_PYR = {
    'F': 14 + 12, 'L': 14 + 13, 'S': 14 + 17, 'Y': 14 + 12,
    'C': 14 + 33, 'W': 14 + 12, 'H': 14 + 12, 'Q': 14 + 14,
    'R': 14 + 14, 'P': 27 + 14, 'stop': 0,
}

# δ and ζ weights from Table 23 (first purine context)
DELTA = {
    'I': 15, 'M': 32, 'T': 0, 'N': 32, 'K': 14, 'S': 0,
    'R': 14, 'V': 0, 'A': 0, 'D': 33, 'E': 12, 'G': 0,
}
ZETA = {
    'I': 0, 'M': 0, 'T': 0, 'N': 0, 'K': 16, 'S': 0,
    'R': 12, 'V': 0, 'A': 0, 'D': 0, 'E': 0, 'G': 0,
}

ARS1 = ['V', 'C', 'L', 'I', 'Q', 'E', 'M', 'R', 'Y', 'W']
ARS2 = ['G', 'A', 'S', 'P', 'T', 'N', 'D', 'K', 'H', 'F']

BASES = ['T', 'C', 'A', 'G']


def side_chain(aa):
    if aa == 'stop':
        return 0
    sc = AMINO_ACIDS[aa][1]
    return 41 if aa == 'P' else sc


def full_weight(aa):
    return 0 if aa == 'stop' else AMINO_ACIDS[aa][2]


def full_weight_stop74(aa):
    return 74 if aa == 'stop' else AMINO_ACIDS[aa][2]


def is_homogeneous(b1, b2):
    products = [GENETIC_CODE[b1 + b2 + b3] for b3 in BASES]
    return len(set(products)) == 1


def sm_counting_products(columns):
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


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")
    return condition


# =============================================================================
# VERIFICATION FUNCTIONS
# =============================================================================

def verify_data_foundation():
    print("\n=== DATA FOUNDATION ===")
    p = True
    p &= check("20 amino acids", len(AMINO_ACIDS) == 20)
    for aa, (name, sc, fw) in AMINO_ACIDS.items():
        cp = fw - sc
        if aa == 'P':
            p &= check("Proline constant part = 73", cp == 73, f"got {cp}")
        elif cp != 74:
            p &= check(f"{name} constant part = 74", False, f"got {cp}")
    p &= check("All non-Pro constant parts = 74",
               all(fw - sc == 74 for aa, (_, sc, fw) in AMINO_ACIDS.items() if aa != 'P'))
    p &= check("Proline native side chain = 42", AMINO_ACIDS['P'][1] == 42)
    p &= check("Proline activated side chain = 41", side_chain('P') == 41)
    return p


def verify_level1_rumer():
    print("\n=== LEVEL 1: RUMER SYMMETRY ===")
    p = True

    octet_I = [(b1, b2) for b1 in BASES for b2 in BASES if is_homogeneous(b1, b2)]
    octet_II = [(b1, b2) for b1 in BASES for b2 in BASES if not is_homogeneous(b1, b2)]

    p &= check("Octet I has 8 homogeneous columns", len(octet_I) == 8,
               f"columns: {[''.join(c) for c in octet_I]}")
    p &= check("Octet II has 8 heterogeneous columns", len(octet_II) == 8)

    # Rumer transformation maps I <-> II
    rumer = {'T': 'G', 'G': 'T', 'C': 'A', 'A': 'C'}
    oi_set, oii_set = set(octet_I), set(octet_II)
    mapped_ok = all((rumer[b1], rumer[b2]) in oii_set for b1, b2 in octet_I)
    mapped_ok &= all((rumer[b1], rumer[b2]) in oi_set for b1, b2 in octet_II)
    p &= check("Rumer transformation maps Octet I <-> Octet II", mapped_ok)

    # 12 unique calligrams
    def make_calligram(perm):
        return frozenset(
            (i, j) for i, b1 in enumerate(perm)
            for j, b2 in enumerate(perm) if is_homogeneous(b1, b2)
        )

    def rotate_180(cal):
        return frozenset((3 - i, 3 - j) for i, j in cal)

    def canonical(cal):
        return tuple(sorted(min(sorted(cal), sorted(rotate_180(cal)))))

    all_perms = list(permutations(BASES))
    calligrams = {perm: make_calligram(perm) for perm in all_perms}
    unique = {canonical(cal) for cal in calligrams.values()}
    p &= check("12 unique calligrams", len(unique) == 12, f"got {len(unique)}")

    # Connected calligrams
    def is_connected_region(cells):
        if not cells:
            return True
        cells = set(cells)
        start = next(iter(cells))
        visited = {start}
        queue = [start]
        while queue:
            ci, cj = queue.pop()
            for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nb = (ci + di, cj + dj)
                if nb in cells and nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
        return len(visited) == len(cells)

    all_cells = {(i, j) for i in range(4) for j in range(4)}
    connected_perms = [
        ''.join(perm) for perm in all_perms
        if is_connected_region(calligrams[perm])
        and is_connected_region(all_cells - calligrams[perm])
    ]
    connected_canonical = {
        canonical(calligrams[tuple(s)]) for s in connected_perms
    }

    p &= check("4 connected calligrams (2 unique under 180° rotation)",
               len(connected_perms) == 4 and len(connected_canonical) == 2,
               f"connected: {connected_perms}")
    p &= check("CTGA is connected", 'CTGA' in connected_perms)
    p &= check("ATGC is connected", 'ATGC' in connected_perms)

    # Group closure under R, R1, R2
    r_full = {'T': 'G', 'G': 'T', 'C': 'A', 'A': 'C'}
    r1 = {'T': 'G', 'G': 'T', 'C': 'C', 'A': 'A'}
    r2 = {'T': 'T', 'G': 'G', 'C': 'A', 'A': 'C'}
    conn_set = {tuple(s) for s in connected_perms}
    group_ok = all(
        tuple(t[b] for b in perm) in conn_set
        for perm in conn_set for t in [r_full, r1, r2]
    )
    p &= check("Connected calligrams closed under R, R1, R2", group_ok)
    return p


def verify_level2_full_weights():
    print("\n=== LEVEL 2: FULL WEIGHT SIGNATURES ===")
    p = True

    octet_I_cols = [(b1, b2) for b1 in BASES for b2 in BASES if is_homogeneous(b1, b2)]
    octet_II_cols = [(b1, b2) for b1 in BASES for b2 in BASES if not is_homogeneous(b1, b2)]

    # Each cell counted separately
    total_I = sum(full_weight(GENETIC_CODE[b1 + b2 + b3])
                  for b1, b2 in octet_I_cols for b3 in BASES)
    total_II = sum(full_weight(GENETIC_CODE[b1 + b2 + b3])
                   for b1, b2 in octet_II_cols for b3 in BASES)
    total_II_74 = sum(full_weight_stop74(GENETIC_CODE[b1 + b2 + b3])
                      for b1, b2 in octet_II_cols for b3 in BASES)

    p &= check("Octet I total weight = 3700", total_I == 3700, f"got {total_I}")
    p &= check("Octet II total weight (stop=0) = 4218", total_II == 4218, f"got {total_II}")
    p &= check("3700 % 74 == 0", 3700 % 74 == 0, f"3700/74 = {3700 // 74}")
    p &= check("4218 % 74 == 0", 4218 % 74 == 0, f"4218/74 = {4218 // 74}")
    p &= check("4218 = 4×999 + 222", 4218 == 4 * 999 + 222)
    p &= check("Octet II total weight (stop=74) = 4440", total_II_74 == 4440, f"got {total_II_74}")
    p &= check("GCD(3700, 4440) = 740", gcd(3700, 4440) == 740)
    p &= check("4440 = 4×999 + 444", 4440 == 4 * 999 + 444)
    p &= check("740 = 2×10×37", 740 == 2 * 10 * 37)
    p &= check("20 = 2×10", 20 == 2 * 10)
    p &= check("74 = 2×37", 74 == 2 * 37)

    # Sorted by weight -> first bases
    octet_I_prods = [(full_weight(GENETIC_CODE[b1 + b2 + 'T']), b1, GENETIC_CODE[b1 + b2 + 'T'])
                     for b1, b2 in octet_I_cols]
    octet_I_prods.sort()
    first_bases = ''.join(x[1] for x in octet_I_prods)
    p &= check("Octet I sorted by weight, first bases = GGTCGACC",
               first_bases == 'GGTCGACC',
               f"got {first_bases} ({[(x[2], x[0]) for x in octet_I_prods]})")

    comp = {'G': 'C', 'C': 'G', 'T': 'A', 'A': 'T'}
    n = len(first_bases)
    mirror = all(first_bases[i] == comp[first_bases[n - 1 - i]] for i in range(n))
    p &= check("GGTCGACC is mirror-complementary", mirror)
    return p


def verify_level3_signature_3_1():
    print("\n=== LEVEL 3, Signature 3.1: Octets side chains ===")
    p = True

    oi_cols = [(b1, b2) for b1 in BASES for b2 in BASES if is_homogeneous(b1, b2)]
    oii_cols = [(b1, b2) for b1 in BASES for b2 in BASES if not is_homogeneous(b1, b2)]

    oi_prods = sm_counting_products(oi_cols)
    sc_I = sum(side_chain(aa) for aa in oi_prods)
    p &= check("Octet I side chains = 333", sc_I == 333, f"got {sc_I}")
    p &= check("333 = 37×9 = 37×3²", 333 == 37 * 9)

    cp_I = len(oi_prods) * 74  # all have cp=74 after proline key
    p &= check("Octet I constant parts = 592", cp_I == 592, f"got {cp_I}")
    p &= check("592 = 37×16 = 37×4²", 592 == 37 * 16)

    fw_I = sum(side_chain(aa) + 74 for aa in oi_prods)
    p &= check("Octet I full weights = 925", fw_I == 925, f"got {fw_I}")
    p &= check("925 = 37×25 = 37×5²", 925 == 37 * 25)
    p &= check("Egyptian triangle: 3² + 4² = 5²", 9 + 16 == 25)

    oii_prods = sm_counting_products(oii_cols)
    sc_II = sum(side_chain(aa) for aa in oii_prods)
    p &= check("Octet II side chains = 1110", sc_II == 1110, f"got {sc_II}")
    p &= check("1110 = 999 + 111", 1110 == 999 + 111)

    cp_II = sum(0 if aa == 'stop' else 74 for aa in oii_prods)
    p &= check("Octet II constant parts (stop=0) = 1110", cp_II == 1110, f"got {cp_II}")
    return p


def verify_level3_signature_3_2():
    print("\n=== LEVEL 3, Signature 3.2: Codons with 2 identical bases ===")
    p = True

    def has_exactly_one_pair(codon):
        b1, b2, b3 = codon
        return sum(a == b for a, b in [(b1, b2), (b1, b3), (b2, b3)]) == 1

    def pair_base(codon):
        b1, b2, b3 = codon
        if b1 == b2: return b1
        if b1 == b3: return b1
        return b2

    codons = [c for c in GENETIC_CODE if has_exactly_one_pair(c)]
    p &= check("36 codons with exactly 2 identical bases", len(codons) == 36)

    pyr = [c for c in codons if pair_base(c) in 'TC']
    pur = [c for c in codons if pair_base(c) in 'AG']
    sc_pyr = sum(side_chain(GENETIC_CODE[c]) for c in pyr)
    sc_pur = sum(side_chain(GENETIC_CODE[c]) for c in pur)
    p &= check("Pyrimidine-pair half: side chains = 999", sc_pyr == 999, f"got {sc_pyr}")
    p &= check("Purine-pair half: side chains = 999", sc_pur == 999, f"got {sc_pur}")
    return p


def verify_level3_signature_3_3():
    print("\n=== LEVEL 3, Signature 3.3: Purine subgroups each = 333 ===")
    p = True

    def has_exactly_one_pair(codon):
        b1, b2, b3 = codon
        return sum(a == b for a, b in [(b1, b2), (b1, b3), (b2, b3)]) == 1

    def pair_base(codon):
        b1, b2, b3 = codon
        if b1 == b2: return b1
        if b1 == b3: return b1
        return b2

    pur_codons = [c for c in GENETIC_CODE if has_exactly_one_pair(c) and pair_base(c) in 'AG']

    def classify(codon):
        b1, b2, b3 = codon
        if b1 == b2:  # adjacent pair at 1,2
            return 1 if b1 == 'A' else 2
        elif b2 == b3:  # adjacent pair at 2,3
            return 1 if b2 == 'A' else 2
        else:  # separated pair at 1,3
            return 3

    groups = {1: [], 2: [], 3: []}
    for c in pur_codons:
        groups[classify(c)].append(c)

    for g in [1, 2, 3]:
        s = sum(side_chain(GENETIC_CODE[c]) for c in groups[g])
        p &= check(f"Purine subgroup {g}: side chains = 333", s == 333, f"got {s}")
    return p


def verify_level3_signature_3_4():
    print("\n=== LEVEL 3, Signature 3.4: Tables 8 & 9 ===")
    p = True

    # Codons from Tables 8 and 9 as explicitly listed in the paper
    table8 = ['TTC', 'CCT', 'AAT', 'AAC', 'GGT', 'GGC',
              'CTT', 'TCC', 'TAA', 'CAA', 'TGG', 'CGG',
              'TCT', 'CTC', 'ATA', 'ACA', 'GTG', 'GCG']
    table9_r1 = ['TTA', 'TTG', 'CCA', 'CCG', 'AAG', 'GGA']
    table9_r2 = ['ATT', 'GTT', 'ACC', 'GCC', 'GAA', 'AGG']
    table9_r3 = ['TAT', 'TGT', 'CAC', 'CGC', 'AGA', 'GAG']
    table9 = table9_r1 + table9_r2 + table9_r3

    s8 = sum(side_chain(GENETIC_CODE[c]) for c in table8)
    s9 = sum(side_chain(GENETIC_CODE[c]) for c in table9)
    r2 = sum(side_chain(GENETIC_CODE[c]) for c in table9_r2)

    p &= check("Table 8 sum = 888", s8 == 888, f"got {s8}")
    p &= check("Table 9 sum = 1110 (999+111)", s9 == 1110, f"got {s9}")
    p &= check("Table 9 row 2 sum = 333", r2 == 333, f"got {r2}")

    # Verify completeness: tables 8+9 = all 36 codons with exactly 2 identical bases
    all_36 = set()
    for codon in GENETIC_CODE:
        b1, b2, b3 = codon
        if sum(a == b for a, b in [(b1, b2), (b1, b3), (b2, b3)]) == 1:
            all_36.add(codon)
    p &= check("Tables 8+9 cover all 36 codons", set(table8 + table9) == all_36)
    return p


def verify_level3_signature_3_5():
    print("\n=== LEVEL 3, Signature 3.5: First pyrimidine/purine total weights ===")
    p = True

    pyr_cols = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in 'TC']
    pur_cols = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in 'AG']

    pyr_total = sum(full_weight_stop74(aa) for aa in sm_counting_products(pyr_cols))
    pur_total = sum(full_weight_stop74(aa) for aa in sm_counting_products(pur_cols))

    p &= check("First pyrimidine total weight (stop=74, SM) = 1776",
               pyr_total == 1776, f"got {pyr_total}")
    p &= check("1776 = 999 + 777", 1776 == 999 + 777)
    p &= check("First purine total weight (SM) = 1517",
               pur_total == 1517, f"got {pur_total}")
    p &= check("1517 = 41 × 37", 1517 == 41 * 37)
    return p


def verify_level3_signature_3_7():
    print("\n=== LEVEL 3, Signature 3.7: Rumer pairs M=(T,G) and K=(C,A) ===")
    p = True

    M, K = 'TG', 'CA'

    def sc_sum(cols):
        return sum(side_chain(aa) for aa in sm_counting_products(cols))

    m1 = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in M]
    k1 = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in K]
    m2 = [(b1, b2) for b1 in BASES for b2 in BASES if b2 in M]
    k2 = [(b1, b2) for b1 in BASES for b2 in BASES if b2 in K]
    mm = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in M and b2 in M]
    kk = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in K and b2 in K]

    p &= check("First base in M: side chains = 654", sc_sum(m1) == 654, f"got {sc_sum(m1)}")
    p &= check("First base in K: side chains = 789", sc_sum(k1) == 789, f"got {sc_sum(k1)}")
    p &= check("Second base in M: side chains = 789", sc_sum(m2) == 789, f"got {sc_sum(m2)}")
    p &= check("Second base in K: side chains = 654", sc_sum(k2) == 654, f"got {sc_sum(k2)}")
    p &= check("Both bases in M: side chains = 369", sc_sum(mm) == 369, f"got {sc_sum(mm)}")
    p &= check("Both bases in K: side chains = 369", sc_sum(kk) == 369, f"got {sc_sum(kk)}")
    return p


def verify_level4_chain_links():
    print("\n=== LEVEL 4: CHAIN LINKS ===")
    p = True

    oi_cols = [(b1, b2) for b1 in BASES for b2 in BASES if is_homogeneous(b1, b2)]
    oii_cols = [(b1, b2) for b1 in BASES for b2 in BASES if not is_homogeneous(b1, b2)]

    # β-level Octet I = 111 (Table 20)
    beta_I = sum(BETA[aa] for aa in sm_counting_products(oi_cols))
    p &= check("β-level Octet I = 111", beta_I == 111, f"got {beta_I}")

    # β+γ levels Octet II = 444 (Table 21)
    bg_II = sum(BETA[aa] + GAMMA_OCTET_II.get(aa, 0)
                for aa in sm_counting_products(oii_cols))
    p &= check("β+γ levels Octet II = 444", bg_II == 444, f"got {bg_II}")

    # First pyrimidine β+γ = 333 (Table 22 — uses its own β+γ values)
    pyr_cols = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in 'TC']
    bg_pyr = sum(BETA_GAMMA_FIRST_PYR.get(aa, BETA[aa])
                 for aa in sm_counting_products(pyr_cols))
    p &= check("First pyrimidine β+γ = 333", bg_pyr == 333, f"got {bg_pyr}")

    # First purine β+δ+ζ = 333 (Table 23)
    pur_cols = [(b1, b2) for b1 in BASES for b2 in BASES if b1 in 'AG']
    bdz_pur = sum(BETA[aa] + DELTA.get(aa, 0) + ZETA.get(aa, 0)
                  for aa in sm_counting_products(pur_cols))
    p &= check("First purine β+δ+ζ = 333", bdz_pur == 333, f"got {bdz_pur}")
    return p


def verify_level5_ars():
    print("\n=== LEVEL 5: ARS SYMMETRIES ===")
    p = True

    # ARS numbering: sorted by side chain weight
    ars1_sorted = sorted(ARS1, key=lambda aa: AMINO_ACIDS[aa][1])
    ars2_sorted = sorted(ARS2, key=lambda aa: AMINO_ACIDS[aa][1])

    ars_num = {}
    for i, aa in enumerate(ars1_sorted):
        ars_num[aa] = i + 1
    for i, aa in enumerate(ars2_sorted):
        ars_num[aa] = -(10 - i)

    # Calligram-B from Table 26: each amino acid assigned to exactly one first-base group
    # The paper assigns each AA to the first base where it appears in a unique/primary column
    calligram_b = {
        'C': ['P', 'L', 'Q', 'H', 'R'],
        'T': ['S', 'C', 'F', 'Y', 'W'],
        'A': ['T', 'I', 'N', 'K', 'M'],
        'G': ['G', 'A', 'V', 'D', 'E'],
    }

    # Verify all 20 amino acids appear exactly once
    all_aa = []
    for row in calligram_b.values():
        all_aa.extend(row)
    p &= check("Calligram-B has all 20 amino acids", sorted(all_aa) == sorted(AMINO_ACIDS.keys()))

    # Verify each row is sorted by side chain weight
    for b1, row in calligram_b.items():
        weights = [AMINO_ACIDS[aa][1] for aa in row]
        p &= check(f"Calligram-B row {b1} sorted by weight",
                   weights == sorted(weights), f"{list(zip(row, weights))}")

    # Column sums
    base_order = ['C', 'T', 'A', 'G']
    col_sums = [0] * 5
    for b1 in base_order:
        for j, aa in enumerate(calligram_b[b1]):
            col_sums[j] += ars_num[aa]

    p &= check("Calligram-B column sums = [-31, 0, 0, 0, +31]",
               col_sums == [-31, 0, 0, 0, 31], f"got {col_sums}")

    print("    Calligram-B:")
    for b1 in base_order:
        entries = [f"{ars_num[aa]:+d},{aa}" for aa in calligram_b[b1]]
        print(f"      {b1}: {entries}")
    print(f"      Σ:  {col_sums}")
    return p


def main():
    print("=" * 60)
    print("Verification of Panov & Filatov (2024)")
    print("Strange Information Structures of the Genetic Code")
    print("=" * 60)

    ap = True
    ap &= verify_data_foundation()
    ap &= verify_level1_rumer()
    ap &= verify_level2_full_weights()
    ap &= verify_level3_signature_3_1()
    ap &= verify_level3_signature_3_2()
    ap &= verify_level3_signature_3_3()
    ap &= verify_level3_signature_3_4()
    ap &= verify_level3_signature_3_5()
    ap &= verify_level3_signature_3_7()
    ap &= verify_level4_chain_links()
    ap &= verify_level5_ars()

    print("\n" + "=" * 60)
    print("ALL CHECKS PASSED" if ap else "SOME CHECKS FAILED")
    print("=" * 60)
    return 0 if ap else 1


if __name__ == '__main__':
    raise SystemExit(main())
