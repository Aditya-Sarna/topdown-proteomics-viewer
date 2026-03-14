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
# Sources: Unimod (unimod.org), PSI-MOD, common top-down proteomics mods
PTM_DATABASE = {
    # Phosphorylation & signalling
    'Phosphorylation':           79.96633,
    'Sulfation':                 79.95681,

    # Acetylation & acylation
    'Acetylation':               42.01057,
    'N-term Acetylation':        42.01057,
    'Propionylation':            56.02621,
    'Butyrylation':              70.04186,
    'Crotonylation':             68.02621,
    'Succinylation':            100.01604,
    'Malonylation':              86.00039,
    'Glutarylation':            114.03169,
    'Formylation':               27.99491,
    '2-Hydroxyisobutyrylation':  86.03678,
    'Benzoylation':             104.02621,

    # Methylation
    'Methylation':               14.01565,
    'Dimethylation':             28.03130,
    'Trimethylation':            42.04695,
    'Monomethyl ester':          14.01565,

    # Ubiquitin / UBL modifications
    'Ubiquitination (GlyGly)':  114.04293,
    'SUMO-1 (Qqtgg)':           396.23429,
    'ISG15 (LRGG)':             340.18558,

    # Oxidation & reduction
    'Oxidation':                 15.99491,
    'Dioxidation':               31.98983,
    'Carbamidomethylation':      57.02146,
    'Propionamide':              71.03711,

    # Deamidation & hydrolysis
    'Deamidation':                0.98402,
    'Glutamine->pyro-Glu':      -17.02655,

    # Pyro-glu (N-terminal)
    'Pyro-glu from Q':          -17.02655,
    'Pyro-glu from E':          -18.01056,
    'Pyro-carbamidomethyl':     -17.02655,

    # Hydroxylation
    'Hydroxylation':             15.99491,
    'Dihydroxylation':           31.98983,

    # Glycosylation
    'HexNAc (O-GlcNAc)':       203.07937,
    'Hex (O-Mannose)':          162.05282,
    'dHex (Fucose)':            146.05791,
    'NeuAc (Sialic acid)':      291.09542,
    'HexNAc+Hex':               365.13219,
    'HexNAc2':                  406.15874,
    'Core1 O-glycan (HexNAc+Hex)': 365.13219,

    # Lipid / GPI
    'Myristoylation':           210.19837,
    'Palmitoylation':           238.22967,

    # Disulfide / Carbamidomethylation variants
    'Dehydration':              -18.01056,

    # Nitrosylation & other reactive-oxygen
    'Nitrosylation':             28.98983,
    'Nitration':                 44.98508,

    # Carbamylation
    'Carbamylation':             43.00581,

    # DiGly remnant / crosslink
    'Diglycine':                114.04293,
}

# PTM preferred target residues
PTM_TARGET_RESIDUES = {
    'Phosphorylation':          'STY',
    'Sulfation':                'Y',
    'Acetylation':              'K',
    'N-term Acetylation':       '',   # N-terminus
    'Propionylation':           'K',
    'Butyrylation':             'K',
    'Crotonylation':            'K',
    'Succinylation':            'K',
    'Malonylation':             'K',
    'Glutarylation':            'K',
    'Formylation':              'K',
    '2-Hydroxyisobutyrylation': 'K',
    'Benzoylation':             'K',
    'Methylation':              'KR',
    'Dimethylation':            'KR',
    'Trimethylation':           'K',
    'Monomethyl ester':         'DE',
    'Ubiquitination (GlyGly)':  'K',
    'SUMO-1 (Qqtgg)':           'K',
    'ISG15 (LRGG)':             'K',
    'Oxidation':                'MW',
    'Dioxidation':              'MW',
    'Carbamidomethylation':     'C',
    'Propionamide':             'C',
    'Deamidation':              'NQ',
    'Glutamine->pyro-Glu':      'Q',
    'Pyro-glu from Q':          'Q',
    'Pyro-glu from E':          'E',
    'Pyro-carbamidomethyl':     'C',
    'Hydroxylation':            'P',
    'Dihydroxylation':          'WF',
    'HexNAc (O-GlcNAc)':       'ST',
    'Hex (O-Mannose)':          'ST',
    'dHex (Fucose)':            'ST',
    'NeuAc (Sialic acid)':      'ST',
    'HexNAc+Hex':               'ST',
    'HexNAc2':                  'ST',
    'Core1 O-glycan (HexNAc+Hex)': 'ST',
    'Myristoylation':           'G',
    'Palmitoylation':           'C',
    'Dehydration':              'ST',
    'Nitrosylation':            'C',
    'Nitration':                'Y',
    'Carbamylation':            'K',
    'Diglycine':                'K',
}

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O',
}
