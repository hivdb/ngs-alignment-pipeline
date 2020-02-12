#! /usr/bin/env python3

import os
import re
import csv
import sys
import pysam
import click
import numpy as np
from collections import defaultdict, Counter, OrderedDict


DIST_INTERPERSON_CUTOFF = 11
DIST_INTRAPERSON_CUTOFF = 4
DIST_INTRASAMPLE_CUTOFF = 3
MIN_FLOAT = sys.float_info.min


def load_patterns(pattern_csv):
    pattern_csv = csv.DictReader(pattern_csv)
    patterns = []
    for row in pattern_csv:
        index = int(row.pop('Index'))
        samplename = row.pop('SampleName')
        count = int(row.pop('Count'))
        pcnt = float(row.pop('Pcnt'))
        pattern = set()
        for pos, na in row.items():
            if na in ('-', '.'):
                continue
            ref = pos[0]
            pos = int(pos[1:])
            pattern.add((pos, na, ref))
        ptid = re.findall(r'^\d+', samplename)[0]
        patterns.append(
            (index, samplename, ptid, pattern, count, pcnt)
        )
    return patterns


def calc_partial_score(dist_array, cutoff, sign=1):
    partials = []
    total_dist = np.log10(len(dist_array) + MIN_FLOAT)
    for dist in range(0, cutoff):
        partials.append((np.log10(
            np.count_nonzero(dist_array == dist) + MIN_FLOAT
        ) - total_dist) * sign)
    partials.append((
        np.log10(np.count_nonzero(dist_array > dist) + MIN_FLOAT) - total_dist
    ) * sign)
    return sum(partials), partials


def calc_score(dist_interperson, dist_intraperson, dist_intrasample):
    s1, p1 = calc_partial_score(dist_interperson, DIST_INTERPERSON_CUTOFF)
    s2, p2 = calc_partial_score(dist_intraperson, DIST_INTRAPERSON_CUTOFF, -1)
    s3, p3 = calc_partial_score(dist_intrasample, DIST_INTRASAMPLE_CUTOFF, -1)
    return (s1 + s2 + s3, *(p1 + p2 + p3))


def na2int(na):
    return {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
        'ins': 4,
        'del': 5
    }[na]


def mut2int(pos, na, ref):
    return np.int64(pos * 36 + na2int(na) * 6 + na2int(ref))


def iter_rows(all_patterns):
    scores = dist2score(all_patterns)
    for row, score in zip(all_patterns, scores):
        yield (*row, *score)
        # np.count_nonzero(dist_interperson == 0),
        # np.count_nonzero(dist_interperson == 1),
        # np.count_nonzero(dist_interperson == 2),
        # np.count_nonzero(dist_interperson == 3),
        # np.count_nonzero(dist_interperson == 4),
        # np.count_nonzero(dist_interperson == 5),
        # np.count_nonzero(dist_interperson >= 6),
        # np.count_nonzero(dist_intraperson == 0),
        # np.count_nonzero(dist_intraperson == 1),
        # np.count_nonzero(dist_intraperson == 2),
        # np.count_nonzero(dist_intraperson == 3),
        # np.count_nonzero(dist_intraperson == 4),
        # np.count_nonzero(dist_intraperson == 5),
        # np.count_nonzero(dist_intraperson >= 6),
        # np.count_nonzero(dist_intrasample == 0),
        # np.count_nonzero(dist_intrasample == 1),
        # np.count_nonzero(dist_intrasample == 2),
        # np.count_nonzero(dist_intrasample == 3),
        # np.count_nonzero(dist_intrasample == 4),
        # np.count_nonzero(dist_intrasample == 5),
        # np.count_nonzero(dist_intrasample >= 6),
        # np.mean(dist_interperson),
        # *np.quantile(dist_interperson, (0, .25, .5, .75, 1)),
        # np.mean(dist_intrasample),
        # *np.quantile(dist_intrasample, (0, .25, .5, .75, 1))


FLAG_INTERPERSON = np.uint8(0x00)
FLAG_INTRAPERSON = np.uint8(0x01)
FLAG_INTRASAMPLE = np.uint8(0x02)


def dist2score(patterns):
    for i, (idx, samplename, ptid, pattern, *_) in enumerate(patterns):
        dist_interperson = []
        dist_intraperson = []
        dist_intrasample = []
        for j, (idx2, samplename2,
                ptid2, pattern2, *_) in enumerate(patterns):
            if idx == idx2:
                continue
            dist = len(pattern ^ pattern2)
            if ptid == ptid2:
                if samplename == samplename2:
                    dist_intrasample.append(dist)
                else:
                    dist_intraperson.append(dist)
            else:
                dist_interperson.append(dist)
        dist_interperson = np.array(dist_interperson)
        dist_intraperson = np.array(dist_intraperson)
        dist_intrasample = np.array(dist_intrasample)
        yield calc_score(
            dist_interperson,
            dist_intraperson,
            dist_intrasample
        )


@click.command()
@click.argument(
    'pattern_csv', type=click.File('r'))
@click.argument(
    'output_csv', type=click.File('w'))
def main(pattern_csv, output_csv):
    patterns = load_patterns(pattern_csv)
    writer = csv.writer(output_csv)
    writer.writerow([
        'SampleName', 'Index', 'PtID', 'Pattern',
        'Count', 'Pcnt',
        'Score',
        *['OutP_eq{}'.format(d) for d in range(DIST_INTERPERSON_CUTOFF)],
        'OutP_gte{}'.format(DIST_INTERPERSON_CUTOFF),
        *['InP_eq{}'.format(d) for d in range(DIST_INTRAPERSON_CUTOFF)],
        'InP_gte{}'.format(DIST_INTRAPERSON_CUTOFF),
        *['InS_eq{}'.format(d) for d in range(DIST_INTRASAMPLE_CUTOFF)],
        'InS_gte{}'.format(DIST_INTRASAMPLE_CUTOFF),
    ])
    # 'OutP_eq0', 'OutP_eq1', 'OutP_eq2',
    # 'OutP_eq3', 'OutP_eq4', 'OutP_eq5', 'OutP_gte6',
    # 'InP_eq0', 'InP_eq1', 'InP_eq2',
    # 'InP_eq3', 'InP_eq4', 'InP_eq5', 'InP_gte6',
    # 'InS_eq0', 'InS_eq1', 'InS_eq2',
    # 'InS_eq3', 'InS_eq4', 'InS_eq5', 'InS_gte6',
    # 'AvgInterpersonDist',
    # 'IP_q0', 'IP_q25', 'IP_q50', 'IP_q75', 'IP_q100',
    # 'AvgIntrasampleDist',
    # 'IS_q0', 'IS_q25', 'IS_q50', 'IS_q75', 'IS_q100',

    for row in iter_rows(patterns):
        (idx, samplename, ptid, pattern,
         count, pcnt, *score) = row
        pattern = ', '.join(['{2}{0}{1}'.format(*m) for m in sorted(pattern)])
        writer.writerow([
            samplename, idx, ptid, pattern, count, pcnt, *score
        ])
    click.echo('DIST_INTERPERSON_CUTOFF: {}'.format(DIST_INTERPERSON_CUTOFF))
    click.echo('DIST_INTRAPERSON_CUTOFF: {}'.format(DIST_INTRAPERSON_CUTOFF))
    click.echo('DIST_INTRASAMPLE_CUTOFF: {}'.format(DIST_INTRASAMPLE_CUTOFF))


if __name__ == '__main__':
    main()
