"""
Python variable aa_info describing amino acid properties. Each property contains a dict with keys as the single-letter
amino acid abbreviations.
    hydrophobic: string, either 'hydrophobic', 'very hydrophobic', 'neutral', or 'hydrophilic'
    name: string, full lowercase names for amino acids
    ph: string indicating ph, either 'neutral', 'acidic' or 'basic'
    polar: string, either 'polar', 'non-polar' or 'charged'
    rchain: string, either 'aliphatic', 'aromatic' or 'other'
    size: string, either 'small', 'medium', 'large', 'large_long chain)' or 'large_benzene)'
"""
aa_info = {'hydrophobicity': {'A': 'hydrophobic', 'C': 'hydrophobic', 'D': 'hydrophillic', 'E': 'hydrophillic',
                           'F': 'very hydrophobic', 'G': 'neutral', 'H': 'neutral', 'I': 'very hydrophobic',
                           'K': 'hydrophillic', 'L': 'very hydrophobic', 'M': 'very hydrophobic', 'N': 'hydrophillic',
                           'P': 'hydrophillic', 'Q': 'neutral', 'R': 'hydrophillic', 'S': 'neutral', 'T': 'neutral',
                           'V': 'very hydrophobic', 'W': 'very hydrophobic', 'Y': 'hydrophobic'},
           'name': {'A': 'alanine', 'C': 'cysteine', 'D': 'aspartic acid', 'E': 'glutamic acid', 'F': 'phenylalanine',
                    'G': 'glycine', 'H': 'histidine', 'I': 'isoleucine', 'K': 'lysine', 'L': 'leucine',
                    'M': 'methionine', 'N': 'asparagine', 'P': 'proline', 'Q': 'glutamine', 'R': 'arginine',
                    'S': 'serine', 'T': 'threonine', 'V': 'valine', 'W': 'tryptophan', 'Y': 'tyrosine'},
           'ph': {'A': 'neutral', 'C': 'neutral', 'D': 'acidic', 'E': 'acidic', 'F': 'neutral', 'G': 'neutral',
                  'H': 'basic', 'I': 'neutral', 'K': 'basic', 'L': 'neutral', 'M': 'neutral', 'N': 'neutral',
                  'P': 'neutral', 'Q': 'neutral', 'R': 'basic', 'S': 'neutral', 'T': 'neutral', 'V': 'neutral',
                  'W': 'neutral', 'Y': 'neutral'},
           'polarity': {'A': 'non-polar', 'C': 'polar', 'D': 'charged', 'E': 'charged', 'F': 'non-polar', 'G': 'non-polar',
                     'H': 'polar', 'I': 'non-polar', 'K': 'charged', 'L': 'non-polar', 'M': 'polar', 'N': 'polar',
                     'P': 'non-polar', 'Q': 'polar', 'R': 'charged', 'S': 'polar', 'T': 'polar', 'V': 'non-polar',
                     'W': 'polar', 'Y': 'polar'},
           'rchain': {'A': 'aliphatic', 'C': 'other', 'D': 'other', 'E': 'other', 'F': 'aromatic', 'G': 'other',
                      'H': 'other', 'I': 'aliphatic', 'K': 'other', 'L': 'aliphatic', 'M': 'other', 'N': 'other',
                      'P': 'other', 'Q': 'other', 'R': 'other', 'S': 'other', 'T': 'other', 'V': 'aliphatic',
                      'W': 'aromatic', 'Y': 'aromatic'},
           'size': {'A': 'small', 'C': 'medium', 'D': 'medium', 'E': 'medium', 'F': 'large-benzene', 'G': 'small',
                    'H': 'large ', 'I': 'medium', 'K': 'large-longchain', 'L': 'medium', 'M': 'medium',
                    'N': 'medium', 'P': 'large ', 'Q': 'large', 'R': 'large-longchain', 'S': 'medium', 'T': 'medium',
                    'V': 'medium', 'W': 'large-benzene', 'Y': 'large-benzene'}
           }

