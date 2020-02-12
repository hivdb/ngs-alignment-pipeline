#! /usr/bin/env python3

from __future__ import print_function

import re
import sys
import csv
from collections import defaultdict, Counter


def main():
    if len(sys.argv) != 3:
        print('Usage: {} <MATRIX_CSV> <DISTSTAT_CSV>'.format(sys.argv[0]),
              file=sys.stderr)
        exit(1)
    matrix_csv, diststat_csv = sys.argv[1:]
    with open(matrix_csv) as matrix_csv, \
            open(diststat_csv, 'w') as diststat_csv:
        reader = csv.reader(matrix_csv)
        writer = csv.writer(diststat_csv)
        headers = next(reader)
        headers.pop(0)  # remove the "##"
        headers = [h.split('_', 1)[0] for h in headers]
        writer.writerow(['Pattern', 'Bin', 'Count'])

        bincounter = defaultdict(Counter)
        for i1, h1 in enumerate(headers):
            row = next(reader)
            for i2, h2 in enumerate(headers):
                if i2 <= i1:
                    continue
                dist = int(float(row[i2]))
                if h1 == h2:
                    bincounter[h1][dist] += 1
                elif re.findall(r'^\d+', h1)[0] != re.findall(r'^\d+', h2)[0]:
                    bincounter['*'][dist] += 1
                    if not (h1.startswith('27911') or h2.startswith('27911')):
                        bincounter['no27911'][dist] += 1
        for ptn in sorted(set(headers)) + ['*', 'no27911']:
            counter = bincounter[ptn]
            for dist, count in sorted(counter.items()):
                writer.writerow([ptn, dist, count])


if __name__ == '__main__':
    main()
