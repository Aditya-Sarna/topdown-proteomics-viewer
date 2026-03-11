# Monoisotopic residue masses (Da) — standard amino acids
AA_MASSES = {
    'A': 71.03711,  'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
    'U': 150.95363, 'O': 237.14773,
}

WATER  = 18.01056
PROTON = 1.007276
NH3    = 17.02655
CO     = 27.99491

# Common PTMs: name -> monoisotopic mass shift (Da)
PTM_DATABASE = {
    'Phosphorylation':           79.96633,
    'Acetylation':               42.01057,
    'Methylation':               14.01565,
    'Dimethylation':             28.03130,
    'Trimethylation':            42.04695,
    'Ubiquitination (GlyGly)':  114.04293,
    'Oxidation':                 15.99491,
    'Deamidation':                0.98402,
    'Carbamidomethylation':      57.02146,
    'Formylation':               27.99491,
    'Hydroxylation':             15.99491,
    'N-term Acetylation':        42.01057,
    'Pyro-glu from Q':          -17.02655,
    'Pyro-glu from E':          -18.01056,
    'Succinylation':            100.01604,
    'Malonylation':              86.00039,
    'Crotonylation':             68.02621,
    'Nitrosylation':             28.98983,
    'Sulfation':                 79.95681,
}

# PTM preferred target residues
PTM_TARGET_RESIDUES = {
    'Phosphorylation':          'STY',
    'Acetylation':              'K',
    'Methylation':              'KR',
    'Dimethylation':            'KR',
    'Trimethylation':           'K',
    'Ubiquitination (GlyGly)':  'K',
    'Oxidation':                'M',
    'Deamidation':              'NQ',
    'Carbamidomethylation':     'C',
    'Formylation':              'K',
    'Hydroxylation':            'P',
    'N-term Acetylation':       '',   # N-terminus
    'Succinylation':            'K',
    'Malonylation':             'K',
    'Crotonylation':            'K',
}

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O',
}
