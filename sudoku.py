"""
For testing.
"""

from StringIO import StringIO
import unittest

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

g_grids = (g_grid_wikipedia, g_grid_hard, g_grid_min_1)


def create_blocks():
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
    return blocks

def get_sets(values, blocks):
    sets = [set() if v else set(range(1, 10)) for v in values]
    for b in blocks:
        taken = set(values[i] for i in b)
        for i in b:
            sets[i] -= taken
    return sets

def solve(values, blocks):
    if all(values):
        return values
    sets = get_sets(values, blocks)
    n, branch = min((len(sets[i]), i) for i, v in enumerate(values) if not v)
    if not n:
        return None
    for value in sets[branch]:
        values[branch] = value
        result = solve(values[:], blocks)
        if result:
            return result
    return None

def solve_faster(values, sets, blocks, inv_blocks):
    if all(values):
        return values
    n, branch = min((len(sets[i]), i) for i, v in enumerate(values) if not v)
    if not n:
        return None
    for value in sets[branch]:
        values[branch] = value
        next_sets = [set(x) for x in sets]
        for b in inv_blocks[branch]:
            for i in b:
                next_sets[i].discard(value)
        result = solve_faster(values[:], next_sets, blocks, inv_blocks)
        if result:
            return result
    return None

def values_to_string(values):
    out = StringIO()
    for i in range(9):
        arr = values[9*i:9*(i+1)]
        print >> out, ' '.join('.' if not x else str(x) for x in arr)
    return out.getvalue().rstrip()


class TestSudoku(unittest.TestCase):

    def test_sudoku_slow(self):
        for grid in g_grids:
            longline = ''.join(grid.split())
            values = [0 if c == '.' else int(c) for c in longline]
            solution = solve(values, create_blocks())
            self.assertTrue(solution is not None)

    def test_sudoku_faster(self):
        blocks = create_blocks()
        inv_blocks = [[] for i in range(81)]
        for b in blocks:
            for i in b:
                inv_blocks[i].append(b)
        for grid in g_grids:
            longline = ''.join(grid.split())
            values = [0 if c == '.' else int(c) for c in longline]
            sets = [set() if v else set(range(1, 10)) for v in values]
            for b in blocks:
                taken = set(values[i] for i in b)
                for i in b:
                    sets[i] -= taken
            solution = solve_faster(values, sets, blocks, inv_blocks)
            self.assertTrue(solution is not None)

    def test_compatibility(self):
        blocks = create_blocks()
        inv_blocks = [[] for i in range(81)]
        for b in blocks:
            for i in b:
                inv_blocks[i].append(b)
        for grid in g_grids:
            longline = ''.join(grid.split())
            values = [0 if c == '.' else int(c) for c in longline]
            sets = [set() if v else set(range(1, 10)) for v in values]
            for b in blocks:
                taken = set(values[i] for i in b)
                for i in b:
                    sets[i] -= taken
            slow_solution = solve(values, blocks)
            faster_solution = solve_faster(values, sets, blocks, inv_blocks)
            self.assertEqual(slow_solution, faster_solution)


if __name__ == '__main__':
    unittest.main()

