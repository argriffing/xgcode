"""
The TikZ project provides a way to embed graphics into LaTeX
documents in a way that makes fonts and line weights look
like the rest of the text.
The R project provides plotting functions with various drivers,
and one of these drivers creates TikZ figures as output.
This is great, except that a bit too much whitespace
is added as a margin around the figure.
This script modifies a TikZ figure generated by R
so that it does not have these excessive margins.

I have since learned that the real way to solve this problem
is by tweaking the margins manually in the R code
as in the following line.
par(mar=c(4.5,4,0,1))
"""

import sys

def processed_lines(line_source):
    """
    Yield EOL terminated lines.
    @param line_source: a source of EOL terminated lines
    """
    # \draw[color=white,opacity=0] (0,0) rectangle (252.94,289.08);
    # \path[clip] (  0.00,  0.00) rectangle (252.94,289.08);
    white_line_suffix = None
    for line in line_source:
        if white_line_suffix is None:
            if all(x in line for x in ('draw', 'color', 'white', 'opacity', 'rectangle')):
                white_line_suffix = line[line.find('rectangle'):]
                continue
        elif line.startswith('\\path[clip]') and line.endswith(white_line_suffix):
            continue
        yield line
    if white_line_suffix is None:
        raise ValueError('Failed to find the opaque white rectangle')

def main():
    for line in processed_lines(sys.stdin):
        sys.stdout.write(line)

if __name__ == '__main__':
    main()
