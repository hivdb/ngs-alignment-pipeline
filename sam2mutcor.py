#! /usr/bin/env python

import os
import csv
import pysam
import click
import numpy as np
from statistics import mean
from multiprocessing import Process, Queue
from itertools import combinations, product
from collections import defaultdict, Counter

from scipy import stats
from skbio import DistanceMatrix
from skbio.tree import nj

import fastareader

REF_CODON_OFFSET = 0

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 30))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 25))
MUTATION_PCNT_CUTOFF = 0.1
PATTERN_PCNT_IDV_CUTOFF = 0.01
PATTERN_PCNT_ACC_CUTOFF = 0.9
MISSING_POSITION_THRESHOLD = 10

NUM_PROCESSES = int(os.environ.get('NTHREADS', 3))
INPUT_QUEUE = Queue(NUM_PROCESSES)
OUTPUT_QUEUE = Queue()

PRODUCER_CAPACITY = 50000

CHUNKSIZE = 500

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


LANL_NAPCNTS = {}
LANL_CONSENSUS = {}
with open('NAPcnts_Dec10_Group_M_Only.csv') as fp:
    reader = csv.DictReader(fp)
    for row in reader:
        LANL_NAPCNTS[(int(row['Pos']), row['Nuc'])] = float(row['Prev']) / 100
        LANL_CONSENSUS[int(row['Pos'])] = row['LANLCons']


def get_na_counts(seq, qua, aligned_pairs, idx, header):

    # pre-filter
    err = ERR_OK
    if len(seq) < LENGTH_CUTOFF:
        err |= ERR_TOO_SHORT
    if mean(qua) < OVERALL_QUALITY_CUTOFF:
        err |= ERR_LOW_QUAL
    if err:
        return (header, None), err

    nas = defaultdict(list)
    prev_refpos = 0
    prev_seqpos0 = 0
    for seqpos0, refpos0 in aligned_pairs:

        if refpos0 is None:
            # insertion
            refpos = prev_refpos
        else:
            refpos = refpos0 + 1
            prev_refpos = refpos

        if seqpos0 is None:
            # deletion
            n, q = '-', qua[prev_seqpos0]
        else:
            n, q = seq[seqpos0], qua[seqpos0]
            prev_seqpos0 = seqpos0

        if refpos == 0:
            continue

        nas[refpos].append((n, q))

    result_nas = []
    for pos, nqs in sorted(nas.items()):
        qs = [q for _, q in nqs]
        meanq = sum(qs) / len(qs)
        if meanq < SITE_QUALITY_CUTOFF:
            continue
        result_nas.append((pos, ''.join(n for n, _ in nqs)))
    if result_nas:
        lastpos, lastna = result_nas[-1]
        # remove insertion at the end of sequence read
        result_nas[-1] = (lastpos, lastna[0])
    return (header, result_nas), err


def reads_consumer():
    while True:
        chunk = INPUT_QUEUE.get()
        out_chunk = []
        for args in chunk:
            out_chunk.append(get_na_counts(*args))
        OUTPUT_QUEUE.put(out_chunk)


def reads_producer(filename, offset):
    with pysam.AlignmentFile(filename, 'rb') as samfile:
        chunk = []
        limit = offset + PRODUCER_CAPACITY
        for idx, read in enumerate(samfile.fetch()):
            if idx < offset:
                continue
            if idx >= limit:
                break
            seq, qua, aligned_pairs, name = (read.query_sequence,
                                             read.query_qualities,
                                             read.get_aligned_pairs(False),
                                             read.query_name)
            if len(chunk) == CHUNKSIZE:
                INPUT_QUEUE.put(chunk)
                chunk = []
            chunk.append((seq, qua, aligned_pairs, idx + offset, name))
        if chunk:
            INPUT_QUEUE.put(chunk)


def pattern2str(pattern, count):
    return 'muts={}|n={}'.format(
        '_'.join('{}{}'.format(*m) for m in pattern),
        count)


def attach_initref_pos(alnprofile):
    result = []
    pos0 = -1
    for p in alnprofile:
        if p != '+':
            pos0 += 1
        result.append((pos0, p))
    return result


def create_producers(sampath, totalreads):
    offset = 0
    processes = []
    while offset < totalreads:
        producer = Process(target=reads_producer, args=(sampath, offset))
        producer.daemon = True
        processes.append(producer)
        producer.start()
        offset += PRODUCER_CAPACITY
    return processes


def create_consumers():
    for _ in range(NUM_PROCESSES):
        consumer = Process(target=reads_consumer)
        consumer.daemon = True
        consumer.start()


def get_patterns(sampath, initrefnas, alnprofile,
                 refbegin, refend, pos_offset):
    with pysam.AlignmentFile(sampath, 'rb') as samfile:
        # set until_eof=False to exclude unmapped reads
        totalreads = samfile.count(until_eof=False)

    producers = create_producers(sampath, totalreads)

    num_finished = 0
    nacounts = Counter()
    patterns = Counter()
    allnacounts = Counter()
    paired_results = defaultdict(list)
    while num_finished < totalreads:
        out_chunk = OUTPUT_QUEUE.get()
        for (header, my_results), my_err in out_chunk:
            num_finished += 1
            paired_results[header].append((my_results, my_err))
    [p.terminate() for p in producers]

    print('Total Reads:', totalreads)
    print('Pairs:', len(paired_results))
    includedreads = 0
    outofrange = 0
    for pairs in paired_results.values():
        if len(pairs) > 2:
            raise RuntimeError('Too many reads in a pair')
        pattern = set()
        allnas = set()
        include = False
        firstpos = [refend] * 2
        lastpos = [refbegin] * 2
        for j, (results, err) in enumerate(pairs):
            if err & ERR_TOO_SHORT or err & ERR_LOW_QUAL:
                firstpos[j] = -1
                lastpos[j] = -1
                continue
            for refpos, nas in results:
                refpos0 = refpos - 1
                initrefpos0, profile = alnprofile[refpos0]
                refna = initrefnas[initrefpos0]
                initrefpos = initrefpos0 + 1
                if initrefpos < refbegin or initrefpos > refend:
                    continue
                include = True
                rel_initrefpos = initrefpos - pos_offset
                # if rel_initrefpos == 120:
                #     print('120', refpos, refna, rel_initrefpos, profile)
                firstpos[j] = min(rel_initrefpos, firstpos[j])
                lastpos[j] = max(rel_initrefpos, lastpos[j])
                if profile == '+':
                    # for tmp_refpos0 in range(refpos0 - 1, -1, -1):
                    #     _, profile = alnprofile[tmp_refpos0]
                    #     if profile != '+':
                    #         refna = refnas[tmp_refpos0]
                    #         break
                    mut = (rel_initrefpos, 'ins', refna)
                else:
                    mut = (rel_initrefpos,
                           nas[0].replace('-', 'del'),
                           refna)
                allnas.add(mut)
                if refna != nas[0]:
                    pattern.add(mut)
        if include:
            includedreads += 1
            patterns[(tuple(sorted(pattern)),
                      *firstpos, *lastpos)] += 1
            nacounts += Counter(pattern)
            allnacounts += Counter(allnas)
            positions = set()
            for pos, _, refna in allnas:
                positions.add((pos, None, refna))
            nacounts += Counter(positions)
            allnacounts += Counter(positions)
        else:
            outofrange += 1
    print('Included Reads:', includedreads)
    print('Out of Range:', outofrange)
    keepmuts = {}
    for (pos, na, refna), count in nacounts.items():
        if na is None:
            continue
        pcnt = count / nacounts[(pos, None, refna)]
        if pcnt > MUTATION_PCNT_CUTOFF:
            keepmuts[(pos, na, refna)] = pcnt
    keepmutset = set(keepmuts.keys())
    final_patterns = Counter()
    for (pattern, *pp), count in patterns.items():
        pattern = tuple(sorted(set(pattern) & keepmutset))
        final_patterns[(pattern, *pp)] += count
    return final_patterns, keepmuts, allnacounts


def logfold(a, b):
    if a == 0:
        return -100
    elif b == 0:
        return 100
    else:
        return np.log10(a) - np.log10(b)


def get_pair_size(firstpos1, firstpos2, lastpos1, lastpos2):
    if firstpos1 == -1:
        return lastpos2 - firstpos2 + 1
    if firstpos2 == -1:
        return lastpos1 - firstpos1 + 1
    positions = sorted([
        (firstpos1, 0), (firstpos2, 0),
        (lastpos1, 1), (lastpos2, 1)
    ])
    directs = tuple(d for _, d in positions)
    if directs == (0, 0, 1, 1):
        return 1 + positions[3][0] - positions[0][0]
    else:  # if directs == (0, 1, 0, 1):
        return 2 + (
            positions[1][0] - positions[0][0] +
            positions[3][0] - positions[2][0])


def filter_patterns(patterns, muts, remove_nas=True):
    collapsed = Counter()
    total = sum(patterns.values())
    idv_cutoff = total * PATTERN_PCNT_IDV_CUTOFF
    acc_cutoff = total * PATTERN_PCNT_ACC_CUTOFF
    for (pattern, firstpos1, firstpos2,
         lastpos1, lastpos2), count in patterns.items():
        if remove_nas and not all(firstpos1 <= m[0] <= lastpos1 or
                                  firstpos2 <= m[0] <= lastpos2
                                  for m in muts):
            continue
        elif not remove_nas:
            # if get_pair_size(firstpos1, firstpos2,
            #                  lastpos1, lastpos2) < minseqlen:
            #     continue
            pattern = set(pattern)
            for pos, na, ref in muts:
                if not (firstpos1 <= pos <= lastpos1 or
                        firstpos2 <= pos <= lastpos2):
                    pattern.add((pos, '.', ref))
            if sum([na == '.' for _, na, _ in pattern]) > \
                    MISSING_POSITION_THRESHOLD:
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
        nfwriter.writerow(['Pos', 'HXB2Cons', 'LANLCons',
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


def export_tree(patterns, treeoutput):
    result_patterns = []
    for idx, (pattern, count, pcnt) in enumerate(patterns):
        result_patterns.append([
            idx, set(pattern),
            '{}_{:.1f}%'.format(
                '-'.join('{2}{0}{1}'.format(*m) for m in pattern),
                pcnt * 100),
            count])

    patterns = result_patterns
    num_patterns = len(patterns)
    if num_patterns < 3:
        print('();')
        return
    dist_matrix = np.zeros((num_patterns, num_patterns), dtype=float)
    patternstrs = [ptnstr for _, _, ptnstr, _ in patterns]
    for (idx1, ptn1, ptnstr1, c1), (idx2, ptn2, ptnstr2, c2) in \
            combinations(patterns, 2):
        distance = len(ptn1 ^ ptn2)  # xor
        dist_matrix[idx1, idx2] = distance
        dist_matrix[idx2, idx1] = distance
    dist_matrix = DistanceMatrix(dist_matrix, patternstrs)
    tree = nj(dist_matrix)
    print(tree.root_at_midpoint())


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
    for (pattern, firstpos1, firstpos2,
         lastpos1, lastpos2), count in patterns.items():
        neg_muts = {(pos, na, r)
                    for pos, na, r in mutset - set(pattern)
                    if firstpos1 <= pos <= lastpos1 or
                    firstpos2 <= pos <= lastpos2}
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
            # 'indelcontext', 'log10fold a+',
            # 'log10fold b+',
            'percent a', 'percent b',
            'percent ab', 'pearson r',
            # 'prevalence a', 'prevalence b'
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
                # pairswithmajorindel[(mut1, mut2)],
                # logfold(
                #     pairs[(mut1, mut2, 0)],
                #     pairs[(mut1, mut2, 1)]
                # ),
                # logfold(
                #     pairs[(mut2, mut1, 0)],
                #     pairs[(mut1, mut2, 1)]
                # ),
                muts[mut1],
                muts[mut2],
                pairs[(mut1, mut2, 1)] / pairs_total,
                calc_pearsonr(
                    pairs[(mut1, mut2, 1)],
                    pairs[(mut1, mut2, 0)],
                    pairs[(mut2, mut1, 0)],
                    pairs[(mut1, mut2, -1)]
                )[0],
                # LANL_NAPCNTS.get(mut1[:2]),
                # LANL_NAPCNTS.get(mut2[:2])
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
    create_consumers()
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
                _, alnprofile = fastareader.load(lastref)
            # refnas = refnas['sequence']
            alnprofile = alnprofile['sequence']
            alnprofile = attach_initref_pos(alnprofile)

            patterns, muts, nacounts = get_patterns(
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
