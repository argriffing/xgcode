"""
Sudoku solver to test textarea presets.

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

g_grid_arto_inkala = """
8........
..36.....
.7..9.2..
.5...7...
....457..
...1...3.
..1....68
..85...1.
.9....4..
""".strip()

def get_form():
    form_objects = [
            Form.MultiLine('sudoku', 'puzzle', g_grid_wikipedia)]
    return form_objects

def get_form_out():
    return FormOut.Report('solution')

def get_presets():
    presets = [
            Form.Preset(
                'a puzzle from the internet',
                {'sudoku' : g_grid_hard}),
            Form.Preset(
                'a 17 clue puzzle',
                {'sudoku' : g_grid_min_1}),
            Form.Preset(
                'tricky by Arto Inkala',
                {'sudoku' : g_grid_arto_inkala}),
            ]
    return presets

def get_response_content(fs):
    longline = ''.join(fs.sudoku.split())
    if len(longline) != 81:
        raise ValueError('expected a 9x9 grid')
    values = [int(c) if c.isdigit() else 0 for c in longline]
    covers = precompute_covers()
    full = set(range(1, 10))
    sets = [full - set(values[i] for i in cover) for cover in covers]
    solution = solve(values, sets, covers)
    if solution is None:
        return 'no solution was found'
    else:
        return values_to_string(solution)

def precompute_covers():
    blocks = []
    for i in range(9):
        blocks.append([i*9 + j for j in range(9)])
        blocks.append([j*9 + i for j in range(9)])
    for i in range(3):
        for j in range(3):
            block = []
            for k in range(3):
                for l in range(3):
                    block.append((i*3+k)*9 + (j*3+l))
            blocks.append(block)
    covers = [set() for i in range(81)]
    for b in blocks:
        for i in b:
            covers[i].update(b)
    return covers

def solve(values, sets, covers):
    """
    @param values: the partial solution
    @param sets: the set of available numbers at each position
    @param covers: the sphere of influence of each position
    """
    if all(values):
        return values
    n, branch = min((len(sets[i]), i) for i, v in enumerate(values) if not v)
    if not n:
        return None
    for value in sets[branch]:
        values[branch] = value
        next_sets = [set(x) for x in sets]
        for i in covers[branch]:
            next_sets[i].discard(value)
        result = solve(values[:], next_sets, covers)
        if result:
            return result
    return None

def values_to_string(values):
    out = StringIO()
    for i in range(9):
        arr = values[9*i:9*(i+1)]
        print >> out, ' '.join('.' if not x else str(x) for x in arr)
    return out.getvalue().rstrip()

