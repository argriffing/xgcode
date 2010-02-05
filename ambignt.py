"""
Deal with ambiguous nucleotides.
"""

g_resolve_nt = {
        'A' : ['A'],
        'C' : ['C'],
        'G' : ['G'],
        'T' : ['T'],
        'R' : ['A', 'G'],
        'Y' : ['C', 'T'],
        'S' : ['G', 'C'],
        'W' : ['A', 'T'],
        'K' : ['G', 'T'],
        'M' : ['A', 'C'],
        'B' : ['C', 'G', 'T'],
        'D' : ['A', 'G', 'T'],
        'H' : ['A', 'C', 'T'],
        'V' : ['A', 'C', 'G'],
        'N' : ['A', 'C', 'G', 'T']}
