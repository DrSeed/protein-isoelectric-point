# Protein net charge, isoelectric point and average molecular weight from sequence.
PKA_SETS = {
    "EMBOSS": {"Nterm": 8.6, "Cterm": 3.6, "K": 10.8, "R": 12.5, "H": 6.5, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1},
    "Lehninger": {"Nterm": 9.69, "Cterm": 2.34, "K": 10.5, "R": 12.4, "H": 6.0, "D": 3.86, "E": 4.25, "C": 8.33, "Y": 10.0},
}
POSITIVE = ("K", "R", "H")
NEGATIVE = ("D", "E", "C", "Y")
RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388, "E": 129.1155,
    "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782, "T": 101.1051, "W": 186.2132,
    "Y": 163.1760, "V": 99.1326,
}
WATER = 18.01528


def net_charge(seq, ph, pka_set="EMBOSS"):
    pka = PKA_SETS[pka_set]
    pos = 1 / (1 + 10 ** (ph - pka["Nterm"]))
    neg = 1 / (1 + 10 ** (pka["Cterm"] - ph))
    for aa in POSITIVE:
        pos += seq.count(aa) / (1 + 10 ** (ph - pka[aa]))
    for aa in NEGATIVE:
        neg += seq.count(aa) / (1 + 10 ** (pka[aa] - ph))
    return pos - neg


def isoelectric_point(seq, pka_set="EMBOSS", tol=1e-4):
    # Net charge falls monotonically with pH, so bisection always finds the zero crossing.
    lo, hi = 0.0, 14.0
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if net_charge(seq, mid, pka_set) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def molecular_weight(seq):
    return sum(RESIDUE_MASS[aa] for aa in seq) + WATER
