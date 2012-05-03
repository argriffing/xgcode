"""
A sudoku solver for testing presets and timeouts.

Solve some sudokus using a plain depth first search
that branches on the grid square with the smallest number
of allowed entries, given the partial solution at each stage.
The default puzzle is the one on the wikipedia sudoku page.
A 17-clue preset is provided
from a collection of minimum sudoku configurations,
and I think it has been proved that no 16-clue puzzle
has a unique solution.
"""

from StringIO import StringIO

import Form
import FormOut

g_grid_wikipedia = """
53..7....
6..195...
.98....6.
8...6...3
4..8.3..1
7...2...6
.6....28.
...419..5
....8..79
""".strip()

g_grid_hard = """
...2...63
3....54.1
..1..398.
.......9.
...538...
.3.......
.263..5..
5.37....8
47...1...
""".strip()

g_grid_min_1 = """
.......1.
4........
.2.......
....5.4.7
..8...3..
..1.9....
3..4..2..
.5.1.....
...8.6...
""".strip()

def get_form():
    form_objects = [
            Form.MultiLine('sudoku',
                'incomplete sudoku grid', g_grid_wikipedia)]
    return form_objects

def get_form_out():
    return FormOut.Report('solution')

def get_presets():
    presets = [
            Form.Preset(
                'arbitrary puzzle labeled hard found using google',
                {'sudoku' : g_grid_hard}),
            Form.Preset(
                'the first 17 clue puzzle from a collection',
                {'sudoku' : g_grid_min_1})
            ]
    return presets

def get_response_content(fs):
    longline = ''.join(fs.sudoku.split())
    values = [0 if c == '.' else int(c) for c in longline]
    if len(values) != 81:
        raise ValueError('expected a 9x9 grid')
    # define the groups of variables associated with the constraints
    blocks_per_index = create_blocks_per_index()
    # define the set of initial choices for each entry
    sets = [set() if v else set(range(1, 10)) for v in values]
    for s, blocks in zip(sets, blocks_per_index):
        s -= set(values[i] for b in blocks for i in b)
    # search for the solution
    solution = solve(values, sets, blocks_per_index)
    if solution is None:
        return 'no solution was found'
    else:
        return values_to_string(solution)

def create_blocks_per_index():
    blocks = []
    # define 1x9 and 9x1 blocks
    for i in range(9):
        blocks.append([i*9 + j for j in range(9)])
        blocks.append([j*9 + i for j in range(9)])
    # define 3x3 blocks
    for i in range(3):
        for j in range(3):
            block = []
            for k in range(3):
                for l in range(3):
                    block.append((i*3+k)*9 + (j*3+l))
            blocks.append(block)
    # define the blocks for each index
    blocks_per_index = [[] for i in range(81)]
    for b in blocks:
        for i in b:
            blocks_per_index[i].append(b)
    return blocks_per_index

def solve(values, sets, blocks_per_index):
    if all(values):
        return values
    # find the variable with the smallest branching factor
    n, branch = min((len(sets[i]), i) for i, v in enumerate(values) if not v)
    if not n:
        return None
    # recursively solve the puzzle
    for value in sets[branch]:
        values[branch] = value
        next_sets = [set(x) for x in sets]
        for b in blocks_per_index[branch]:
            for i in b:
                next_sets[i].discard(value)
        result = solve(values[:], next_sets, blocks_per_index)
        if result:
            return result
    return None

def values_to_string(values):
    out = StringIO()
    for i in range(9):
        arr = values[9*i:9*(i+1)]
        print >> out, ' '.join('.' if not x else str(x) for x in arr)
    return out.getvalue().rstrip()

