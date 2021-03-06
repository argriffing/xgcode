"""
"""


def get_folded_hamming_distance(foo):
    return foo

# toggle each row that starts with a one
snps = [
        [0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0],
        [0,0,1,0,1,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,1,1,0,1,1,1,0,0,1,1,1,0,0,1,0,1,1,1],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1],
        [0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,1,0,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,1,1,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [1,1,0,1,0,0,1,0,1,0,1,1,0,0,0,1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,0,1],
        [0,0,1,0,1,1,1,1,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,1,1,0,1,0,0,1,0,1,0,0,0],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,1,1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,1,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0],
        [0,0,1,1,1,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,1],
        [0,0,1,0,1,1,1,1,1,0,0,0,1,1,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,1,1,1],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0],
        [0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,1,1,1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,1],
        [0,0,1,1,1,1,0,1,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,1,1,0,1,0,1,1,0]]


