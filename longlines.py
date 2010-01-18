"""
Detect lines that are longer than the recommended length.
"""

import os

def main():
    for filename in os.listdir('.'):
        if filename.endswith('.py'):
            f = open(filename, 'r')
            for line_number, line in enumerate(f):
                if len(line) > 75:
                    print filename, line_number + 1, line.rstrip()
            f.close()

if __name__ == '__main__':
    main()
