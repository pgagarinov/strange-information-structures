# The Rumer Octets Are Weight-Sorted

**A previously unnoticed property of the genetic code's information structures**

## Background

The universal genetic code table has 16 columns, defined by the first two bases of each codon. Yuri Rumer discovered that these columns split into two groups of 8: **Octet I** (homogeneous columns, where all four codons encode the same amino acid) and **Octet II** (heterogeneous columns, encoding two or more products). This partition is purely structural — it depends only on the degeneracy pattern of the code, not on any chemical property of the amino acids.

Separately, Shcherbak and Makukov (2013), and later Panov and Filatov (2024), showed that summing the side chain weights of amino acids in each octet (using the SM counting rule: each distinct product counted once per column, with the proline activation key 42→41 and stop=0) yields striking results:

- Octet I: **333** = 37 × 9
- Octet II: **1110** = 999 + 111

Both are "SM-numbers" — numbers of the form *m* × 999 + *n* × 111 — which play a central role in the information signatures of the genetic code.

## The observation

When the 16 columns are sorted by their SM-counted side chain weight, the Rumer boundary falls exactly at position 8:

| Rank | Column | Weight | Amino acid(s) | Octet |
|------|--------|--------|---------------|-------|
| 1 | GG | 1 | Gly | I |
| 2 | GC | 15 | Ala | I |
| 3 | TC | 31 | Ser | I |
| 4 | CC | 41 | Pro | I |
| 5 | GT | 43 | Val | I |
| 6 | AC | 45 | Thr | I |
| 7 | CT | 57 | Leu | I |
| 8 | CG | 100 | Arg | I |
| 9 | TA | 107 | Tyr, stop | II |
| 10 | AA | 130 | Asn, Lys | II |
| 11 | AG | 131 | Ser, Arg | II |
| 12 | AT | 132 | Ile, Met | II |
| 13 | GA | 132 | Asp, Glu | II |
| 14 | TT | 148 | Phe, Leu | II |
| 15 | CA | 153 | His, Gln | II |
| 16 | TG | 177 | Cys, stop, Trp | II |

**Every Octet I column is lighter than every Octet II column.** The separation is perfect — zero overlap, zero exceptions. The heaviest Octet I column (CG, Arg, weight 100) is still lighter than the lightest Octet II column (TA, Tyr/stop, weight 107).

## Why this is strange

The Rumer partition is defined by a **logical** property: whether a column's four codons all encode the same product. The weight ranking depends on the **chemistry** of amino acid side chains — the number and type of atoms in each molecule. These are, a priori, unrelated.

Yet the correlation is not merely statistical — it is total. The genetic code assigned its eight simplest amino acids (by side chain mass) to the eight homogeneous columns, and its more complex products to the heterogeneous columns, without a single exception.

This means:

1. **333 is the minimum.** No other choice of 8 columns from the 16 produces a smaller sum. The Rumer octets are pinned to the extremes of the weight distribution.

2. **1110 is the maximum.** Likewise for Octet II.

3. **The total is an SM-number.** 333 + 1110 = 1443 = 999 + 444 = 111 × 13. This was not noted by Panov and Filatov.

## Consequence for probability estimates

Panov and Filatov estimate the probability of both octet sums being divisible by 74 as approximately 1/5476, treating the sums as random integers. But the perfect weight-sorting means the octet sums are not arbitrary — they are the minimum and maximum possible 8-column sums. This concentrates the values at the tails of the distribution rather than the middle, which changes the probability calculation. The octets are structurally constrained to be extreme, making the SM-number property of those extremes a sharper coincidence than a uniform random model would suggest.

## An analogy

Imagine 16 cards, each bearing a number. You partition them into two stacks — not by looking at the numbers, but by a property printed on the back (say, a circle or a square). You then discover that the circle stack contains exactly the 8 smallest numbers and the square stack the 8 largest, with no overlap. And both sums happen to be multiples of 111.

In the genetic code, the "back of the card" is the degeneracy structure (homogeneous vs. heterogeneous), and the "number on the face" is the amino acid side chain mass. Two apparently independent classification systems — one logical, one chemical — produce the same partition.

## References

- Panov, A.D. and Filatov, F.P. (2024). Are the Strange Information Structures of the Genetic Code an Accident or an Artifact? *Evolution: Environmental, Demographic, and Political Risks*, pp. 23–68.
- Shcherbak, V. and Makukov, M. (2013). The "Wow! Signal" of the Terrestrial Genetic Code. *Icarus*, 224, pp. 228–242.
- Rumer, Yu. (2013 [1966]). Systematization of Codons in the Genetic Code. *DAN SSSR*, 183, pp. 222–226.
