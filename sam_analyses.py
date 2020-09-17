#! /usr/bin/env python

import os
import csv
import click
import numpy as np
from itertools import combinations, product
from collections import defaultdict, Counter

from patternutils import find_patterns, filter_patterns

from scipy import stats
from skbio import DistanceMatrix
from skbio.tree import nj

import fastareader

MUTATION_PCNT_CUTOFF = 0.1

LANL_NAPCNTS = {}
LANL_CONSENSUS = {}
with open('NAPcnts_Dec10_Group_M_Only.csv') as fp:
    reader = csv.DictReader(fp)
    for row in reader:
        LANL_NAPCNTS[(int(row['Pos']), row['Nuc'])] = float(row['Prev']) / 100
        LANL_CONSENSUS[int(row['Pos'])] = row['LANLCons']


def export_pattern_table(patterns, muts, patternoutput):
    mutlist = sorted(muts.keys())
    with open(patternoutput, 'w') as patternoutput:
        ptnwriter = csv.writer(patternoutput)
        ptnwriter.writerow(
            ['{2}{0}{1}'.format(*m) for m in mutlist] +
            ['Count', 'Pcnt']
        )
        rows = []
        for pattern, count, pcnt in patterns:
            row = [
                'Y' if m in pattern else ''
                for m in mutlist
            ]
            rows.append([*row, count, pcnt])
        ptnwriter.writerows(sorted(rows))


def export_nucfreq_table(nacounts, nucfreqoutput, preserve_ins_detail):
    with open(nucfreqoutput, 'w') as nucfreqoutput:
        totals = [(m, c) for m, c in nacounts.items()
                  if m[1] is None]
        nfwriter = csv.writer(nucfreqoutput)
        nfwriter.writerow(['Pos', 'HXB2', 'LANLCons',
                           'Nuc', 'Count', 'Total', 'Pcnt', 'Prev'])
        for (pos, _, ref), total in sorted(totals):
            options = ['A', 'C', 'G', 'T']
            if preserve_ins_detail:
                for _pos, na, _ in nacounts.keys():
                    if pos == _pos and na and \
                            len(na) > 1 and na not in ('ins', 'del'):
                        options.append(na)
            options.append('ins')
            options.append('del')
            ins_count = 0
            for na in options:
                count = nacounts[(pos, na, ref)]
                if len(na) > 1 and na not in ('ins', 'del'):
                    ins_count += count
                if na == 'ins':
                    count = ins_count
                nfwriter.writerow([
                    pos, ref,
                    LANL_CONSENSUS.get(pos),
                    na, count, total,
                    count / total,
                    LANL_NAPCNTS.get((pos, na)),
                ])


def calc_pearsonr(both, only_a, only_b, neither):
    total = sum([both, only_a, only_b, neither])
    if total < 2:
        return (None, None)
    x, y = np.zeros(total, int), np.zeros(total, int)
    x[:both] = y[:both] = 1
    x[both:both + only_a] = 1
    y[both + only_a:total - neither] = 1
    return stats.pearsonr(x, y)


def export_mutcor_table(patterns, muts, samplename, mutcoroutput):
    pairs = Counter()
    pairswithmajorindel = Counter()
    mutset = set(muts.keys())
    for (pattern, coverage), count in patterns.items():
        coverage = set(coverage)
        neg_muts = {(pos, na, r)
                    for pos, na, r in mutset - set(pattern)
                    if pos in coverage}
        indelpos = sorted({
            pos for pos, na, r in pattern
            if na in ('ins', 'del')})
        for mut12 in combinations(neg_muts, 2):
            mut1, mut2 = sorted(mut12)
            pairs[(mut1, mut2, -1)] += count
        for mut1, mut2 in product(pattern, neg_muts):
            pairs[(mut1, mut2, 0)] += count
            if any(p <= mut1[0] for p in indelpos):
                pairswithmajorindel[(mut1, mut2)] += count
        for mut12 in combinations(pattern, 2):
            mut1, mut2 = sorted(mut12)
            pairs[(mut1, mut2, 1)] += count
            if any(p <= mut2[0] for p in indelpos):
                pairswithmajorindel[(mut1, mut2)] += count
    with open(mutcoroutput, 'w') as mutcoroutput:
        writer = csv.writer(mutcoroutput)
        writer.writerow([
            'sample name',
            'mutation a', 'mutation b',
            'a+b+', 'a+b-', 'a-b+', 'a-b-',
            'percent a', 'percent b',
            'percent ab', 'pearson r',
        ])
        for mut12 in combinations(mutset, 2):
            mut1, mut2 = sorted(mut12)
            if mut1[0] == mut2[0]:
                continue
            pairs_total = (
                pairs[(mut1, mut2, 1)] +
                pairs[(mut1, mut2, 0)] +
                pairs[(mut2, mut1, 0)] +
                pairs[(mut1, mut2, -1)]
            )
            if pairs_total == 0:
                continue
            writer.writerow([
                samplename,
                '{2}{0}{1}'.format(*mut1),
                '{2}{0}{1}'.format(*mut2),
                pairs[(mut1, mut2, 1)],
                pairs[(mut1, mut2, 0)],
                pairs[(mut2, mut1, 0)],
                pairs[(mut1, mut2, -1)],
                muts[mut1],
                muts[mut2],
                pairs[(mut1, mut2, 1)] / pairs_total,
                calc_pearsonr(
                    pairs[(mut1, mut2, 1)],
                    pairs[(mut1, mut2, 0)],
                    pairs[(mut2, mut1, 0)],
                    pairs[(mut1, mut2, -1)]
                )[0],
            ])


def build_poslookup(mutlist):
    poslookup = defaultdict(list)
    for pos, na, ref in mutlist:
        poslookup[(pos, ref)].append(na)
    return poslookup


def export_all_pattern_table(all_patterns, all_muts, output):
    with open(output, 'w') as output:
        writer = csv.writer(output)
        all_poslookup = build_poslookup(all_muts)
        allpos = sorted(all_poslookup)
        writer.writerow([
            'SampleName', 'Index',
            *['{1}{0}'.format(*p) for p in allpos],
            'Count', 'Pcnt'])
        for idx, (samplename,
                  pattern, count, pcnt) in enumerate(all_patterns):
            poslookup = build_poslookup(pattern)
            row = [samplename, idx + 1]
            row += [
                ''.join(poslookup.get(pos, '-'))
                for pos in allpos
            ]
            writer.writerow([*row, count, pcnt])


def export_tree_for_all(all_patterns, matrixoutput, treeoutput):
    result_patterns = []
    for idx, (samplename, pattern, count, pcnt) in enumerate(all_patterns):
        removes = set()
        for pos, na, ref in pattern:
            if na == '.':
                removes |= {
                    (pos, ntmp, ref)
                    for ntmp in ['A', 'C', 'G', 'T', 'ins', 'del', '.']
                }
        result_patterns.append([
            idx, set(pattern) - removes, removes,
            '{}_{}_{:.1f}%'.format(
                samplename,
                idx + 1,
                pcnt * 100),
            count])

    patterns = result_patterns
    num_patterns = len(patterns)
    if num_patterns < 3:
        with open(treeoutput, 'w') as fp:
            fp.write('();')
        return
    dist_matrix = np.zeros((num_patterns, num_patterns), dtype=float)
    patternstrs = [ptnstr for _, _, _, ptnstr, _ in patterns]
    for (idx1, ptn1, rm1, ptnstr1, c1), (idx2, ptn2, rm2, ptnstr2, c2) in \
            combinations(patterns, 2):
        distance = len((ptn1 - rm2) ^ (ptn2 - rm1))  # xor
        dist_matrix[idx1, idx2] = distance
        dist_matrix[idx2, idx1] = distance
    with open(matrixoutput, 'w') as fp:
        writer = csv.writer(fp)
        writer.writerow(['##', *patternstrs])
        writer.writerows(dist_matrix)
    if True or num_patterns > 10000:
        # TODO: add a switch to this
        # Too many patterns, unable to calculate dist_matrix
        return
    dist_matrix = DistanceMatrix(dist_matrix, patternstrs)
    tree = nj(dist_matrix)
    with open(treeoutput, 'w') as fp:
        fp.write(str(tree.root_at_midpoint()))


@click.command()
@click.argument(
    'initref', type=click.File('r'))
@click.argument(
    'samfolder', type=click.Path(dir_okay=True, exists=True))
@click.option('--ref-range', default='0-99999')
@click.option('--pos-offset', type=int, default=0)
@click.option('--idv-cutoff', type=float, required=True, default=0.01)
@click.option('--acc-cutoff', type=float, required=True, default=0.9)
@click.option('--miss-pos-cutoff', type=int, required=True, default=10)
@click.option(
    '--preserve-ins-detail', is_flag=True,
    help='Preserve insertion detail or just use "ins"')
def main(samfolder, initref, ref_range, pos_offset,
         idv_cutoff, acc_cutoff, miss_pos_cutoff,
         preserve_ins_detail):
    refbegin, refend = [int(r) for r in ref_range.split('-', 1)]
    initrefnas, = fastareader.load(initref)
    initrefnas = initrefnas['sequence']
    all_patterns = []
    all_muts = set()
    for basedir, _, filenames in os.walk(samfolder):
        for sampath in filenames:
            if sampath[-4:].lower() not in ('.sam', '.bam'):
                continue
            sampath = os.path.join(basedir, sampath)
            lastref = os.path.splitext(
                os.path.splitext(sampath)[0]
            )[0] + '.lastref.fas'
            mutcoroutput = os.path.splitext(sampath)[0] + '.mutcor.csv'
            patternoutput = os.path.splitext(sampath)[0] + '.pattern.csv'
            nucfreqoutput = os.path.splitext(sampath)[0] + '.nucfreq.csv'
            with open(lastref) as lastref:
                _, lastrefprofile = fastareader.load(lastref)
            lastrefprofile = lastrefprofile['sequence']

            patterns, muts, nacounts = find_patterns(
                sampath, initrefnas, lastrefprofile,
                refbegin, refend, pos_offset, MUTATION_PCNT_CUTOFF,
                preserve_ins_detail
            )
            filtered_patterns = filter_patterns(
                patterns, muts.keys(),
                idv_cutoff, acc_cutoff, miss_pos_cutoff,
                remove_partials=False)

            export_pattern_table(filtered_patterns, muts, patternoutput)
            export_nucfreq_table(nacounts, nucfreqoutput, preserve_ins_detail)

            samplename = sampath.rsplit('/', 1)[-1].split('_', 1)[0]
            export_mutcor_table(patterns, muts, samplename, mutcoroutput)

            all_patterns.append((samplename, patterns))
            all_muts |= set(muts.keys())
            click.echo('{} processed'.format(sampath))

    all_patterns = [
        (samplename, *ptnresult)
        for samplename, ptns in all_patterns
        for ptnresult in filter_patterns(
            ptns, all_muts,
            idv_cutoff, acc_cutoff, miss_pos_cutoff,
            remove_partials=False
        )
    ]
    export_all_pattern_table(
        all_patterns, all_muts,
        os.path.join(samfolder, 'all_patterns.csv'))
    export_tree_for_all(
        all_patterns,
        os.path.join(samfolder, 'all_patterns.matrix.csv'),
        os.path.join(samfolder, 'all_patterns.tree'))


if __name__ == '__main__':
    main()
