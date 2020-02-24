#! /usr/bin/env python

import os
import re
import csv
import ssw
import click
import numpy as np
from itertools import combinations, product
from collections import defaultdict, Counter

from samreader import iter_paired_reads

from scipy import stats
from skbio import DistanceMatrix
from skbio.tree import nj

import fastareader

MUTATION_PCNT_CUTOFF = 0.1
PATTERN_PCNT_IDV_CUTOFF = 0.01
PATTERN_PCNT_ACC_CUTOFF = 0.9
MISSING_POSITION_THRESHOLD = 10

LANL_NAPCNTS = {}
LANL_CONSENSUS = {}
with open('NAPcnts_Dec10_Group_M_Only.csv') as fp:
    reader = csv.DictReader(fp)
    for row in reader:
        LANL_NAPCNTS[(int(row['Pos']), row['Nuc'])] = float(row['Prev']) / 100
        LANL_CONSENSUS[int(row['Pos'])] = row['LANLCons']


def attach_initref_pos(alnprofile):
    result = []
    pos0 = -1
    for p in alnprofile:
        if p != '+':
            pos0 += 1
        result.append((pos0, p))
    return result


CIGAR_PATTERN = re.compile(r'(\d+)([MNDI])')


def swalign(reference, query):
    aligner = ssw.Aligner()
    alignment = aligner.align(reference=reference, query=query)
    cigar = alignment.cigar
    alignment_profile = []
    offset = 0
    for size, flag in CIGAR_PATTERN.findall(cigar):
        size = int(size)
        if flag in 'MND':
            # for M(atch), N(ull) and D(el) flags
            alignment_profile += list(range(offset, offset + size))
            offset += size
        else:
            # for I(ns) flags
            alignment_profile += [offset] * size
    return alignment_profile


def map_allnas_to_initref(allnas, initrefnas, alnprofile):
    """
    Re-map an iteratively aligned reads to the original reference

    The position numbers of original reference (initref) was lost during
    the iterative alignment. This generator function restores the position
    numbers by using information from alignment profile list (alnprofile).

    This generator produces a tuple of three elements:
      - `position`: The original reference position minus `pos_offset`.
      - `nas`: Nucleic acid notations of this position. The insertion gap is
        represented in higher numbers (>1) of `nas`.
      - `refna`: Nucleic acid notation of the orig. reference at this positon.
    """
    initrefpos_map = defaultdict(list)
    for refpos, nas, _ in allnas:
        refpos0 = refpos - 1
        initrefpos0 = alnprofile[refpos0]
        initrefpos_map[initrefpos0].append(nas)
    for initrefpos0, nas in sorted(initrefpos_map.items()):
        nas = ''.join(nas)
        if len(nas) > 1:
            # remove next-to-insertion deletions
            nas = nas.replace('-', '')
        initrefpos = initrefpos0 + 1
        refna = initrefnas[initrefpos0]
        yield initrefpos, nas, refna


def trim_and_offset_allnas(allnas, refbegin, refend, ref_offset):
    """
    Trim input `allnas` and offset the position by `ref_offset
    """
    for refpos, nas, refna in allnas:
        if refpos < refbegin or refpos > refend:
            continue
        rel_refpos = refpos - ref_offset
        yield rel_refpos, nas, refna


def replace_indel_notations(allnas):
    """
    Replace the indel notations from long NAs and "-" to "ins" and "del"
    """
    for refpos, nas, refna in allnas:
        if len(nas) > 1:
            yield refpos, 'ins', refna
        else:
            yield refpos, nas.replace('-', 'del'), refna


def extract_mutations(allnas):
    """
    Extract mutations (diff from refna)

    Here we assume that refna doesn't have indels
    """
    return ((p, n, r) for p, n, r in allnas if n != r)


def extract_positions(allnas):
    """Extract positions"""
    return (p for p, _, _ in allnas)


def find_patterns(sampath, initrefnas, alnprofile,
                  refbegin, refend, pos_offset):
    patterns = Counter()
    mutcounts = Counter()
    nacounts = Counter()
    for header, allnas in iter_paired_reads(sampath):

        # pre-process allnas
        allnas = map_allnas_to_initref(allnas, initrefnas, alnprofile)
        allnas = trim_and_offset_allnas(allnas, refbegin, refend, pos_offset)
        allnas = set(replace_indel_notations(allnas))

        # extract mutations & positions
        mutations = set(extract_mutations(allnas))
        coverage = set(extract_positions(allnas))

        if allnas:
            # pattern counter ++
            patterns[(tuple(sorted(mutations)),
                      tuple(sorted(coverage)))] += 1

            # mutation counter ++
            mutcounts += Counter(mutations)

            # nucleic acid counter ++
            nacounts += Counter(allnas)

            # position counter ++
            nacounts += Counter({(pos, None, refna)
                                 for pos, _, refna in allnas})

    # find out mutations about MUTATION_PCNT_CUTOFF
    keepmuts = {}
    for (pos, na, refna), count in mutcounts.items():
        pcnt = count / nacounts[(pos, None, refna)]
        if pcnt >= MUTATION_PCNT_CUTOFF:
            keepmuts[(pos, na, refna)] = pcnt
    keepmutset = set(keepmuts.keys())

    # remove mutations below MUTATION_PCNT_CUTOFF from patterns and re-count
    final_patterns = Counter()
    for (pattern, coverage), count in patterns.items():
        pattern = tuple(sorted(set(pattern) & keepmutset))
        final_patterns[(pattern, coverage)] += count
    return final_patterns, keepmuts, nacounts


def filter_patterns(patterns, muts, remove_partials=True):
    """Apply filter rules to patterns generated by `find_patterns`

    Rules:

      - MISSING_POSITION_THRESHOLD: Without remove_partials=True, the
        function allows to keep a partial pattern, if its position coverage
        is above this threshold. The pattern will be furtherly judged by
        the two below rules.
      - PATTERN_PCNT_ACC_CUTOFF: Sorted patterns by percent decendingly,
        then retain the top PATTERN_PCNT_ACC_CUTOFF (0-1.0) of patterns.
      - PATTERN_PCNT_IDV_CUTOFF: Patterns above this threshold are all
        kept no matter to PATTERN_PCNT_ACC_CUTOFF set.

    """
    collapsed = Counter()
    total = sum(patterns.values())
    idv_cutoff = total * PATTERN_PCNT_IDV_CUTOFF
    acc_cutoff = total * PATTERN_PCNT_ACC_CUTOFF
    mutpositions = {p for p, _, _ in muts}
    for (pattern, coverage), count in patterns.items():
        coverage = set(coverage)
        if remove_partials:
            # Partial covered pattern is not allowed
            # if specified remove_partials=True.
            if not coverage.issupperset(mutpositions):
                continue
        else:
            # Partial covered pattern is allowed.
            # Populate missed positions into the pattern
            pattern = set(pattern)
            missed = set()
            for pos, _, ref in muts:
                if pos not in coverage:
                    pattern.add((pos, '.', ref))
                    missed.add(pos)

            # check and skip if the number of missed
            # positions exceeds threshold
            if len(missed) > MISSING_POSITION_THRESHOLD:
                continue
            pattern = tuple(sorted(pattern))

        collapsed[pattern] += count

    results = []
    countsum = 0
    for pattern, count in collapsed.most_common():
        if countsum > acc_cutoff and count < idv_cutoff:
            break
        results.append([pattern, count, count / total])
        countsum += count
    return results


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


def export_nucfreq_table(nacounts, nucfreqoutput):
    with open(nucfreqoutput, 'w') as nucfreqoutput:
        totals = [(m, c) for m, c in nacounts.items()
                  if m[1] is None]
        nfwriter = csv.writer(nucfreqoutput)
        nfwriter.writerow(['Pos', 'HXB2', 'LANLCons',
                           'Nuc', 'Count', 'Total', 'Pcnt', 'Prev'])
        for (pos, _, ref), total in sorted(totals):
            for na in ('A', 'C', 'G', 'T', 'ins', 'del'):
                count = nacounts[(pos, na, ref)]
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
    if PATTERN_PCNT_ACC_CUTOFF == 1:
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
@click.option('--acc-cutoff', type=float, required=True)
def main(samfolder, initref, ref_range, pos_offset, acc_cutoff):
    global PATTERN_PCNT_ACC_CUTOFF
    PATTERN_PCNT_ACC_CUTOFF = acc_cutoff
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
                lastref, _ = fastareader.load(lastref)
            lastref = lastref['sequence']

            # Realign lastref using initrefnas as the reference.
            # The new alignment will be used to restore the initrefnas
            # position numbers.
            alnprofile = swalign(initrefnas, lastref)

            patterns, muts, nacounts = find_patterns(
                sampath, initrefnas, alnprofile,
                refbegin, refend, pos_offset)
            filtered_patterns = filter_patterns(patterns, muts.keys(), False)

            export_pattern_table(filtered_patterns, muts, patternoutput)
            export_nucfreq_table(nacounts, nucfreqoutput)

            samplename = sampath.rsplit('/', 1)[-1].split('_', 1)[0]
            export_mutcor_table(patterns, muts, samplename, mutcoroutput)

            all_patterns.append((samplename, patterns))
            all_muts |= set(muts.keys())
            click.echo('{} processed'.format(sampath))

    all_patterns = [
        (samplename, *ptnresult)
        for samplename, ptns in all_patterns
        for ptnresult in filter_patterns(ptns, all_muts, False)
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
