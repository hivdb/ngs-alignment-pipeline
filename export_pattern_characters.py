#! /usr/bin/env python

import os
import re
import csv
import sys
import pysam
import click
import numpy as np
from collections import defaultdict, Counter, OrderedDict

REF_CODON_OFFSET = 0

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 30))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 25))
MUTATION_PCNT_CUTOFF = 0.1

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


def get_na_counts(read):
    seq, qua, aligned_pairs = (read.query_sequence,
                               read.query_qualities,
                               read.get_aligned_pairs(False))

    # pre-filter
    err = ERR_OK
    if len(seq) < LENGTH_CUTOFF:
        err |= ERR_TOO_SHORT
    if np.mean(qua) < OVERALL_QUALITY_CUTOFF:
        err |= ERR_LOW_QUAL
    if err:
        return None, err

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
    return result_nas, err


def iter_paired_reads(samfile):
    paired_reads = OrderedDict()
    with pysam.AlignmentFile(samfile, 'rb') as samfile:
        for idx, read in enumerate(samfile.fetch()):
            name = read.query_name
            paired_reads.setdefault(name, []).append(read)
    yield from paired_reads.values()


def attach_initref_pos(alnprofile):
    result = []
    pos0 = -1
    for p in alnprofile:
        if p != '+':
            pos0 += 1
        result.append((pos0, p))
    return result


def find_patterns(sampath, initrefnas, alnprofile,
                  refbegin, refend, pos_offset):

    # find out full mutation list for each reads and count by mutations
    nacounts = Counter()
    read_patterns = []
    for pairs in iter_paired_reads(sampath):
        if len(pairs) > 2:
            raise RuntimeError('Too many reads in a pair')
        pattern = set()
        include = False
        firstpos = [refend] * 2
        lastpos = [refbegin] * 2
        for j, read in enumerate(pairs):
            results, err = get_na_counts(read)
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
                # nacounts[(rel_initrefpos, None)] += 1
                if profile == '+':
                    # for tmp_refpos0 in range(refpos0 - 1, -1, -1):
                    #     _, profile = alnprofile[tmp_refpos0]
                    #     if profile != '+':
                    #         refna = refnas[tmp_refpos0]
                    #         break
                    mut = (rel_initrefpos, 'ins', refna)
                elif refna != nas[0]:
                    mut = (rel_initrefpos,
                           nas[0].replace('-', 'del'),
                           refna)
                else:
                    continue
                pattern.add(mut)
                # nacounts[mut] += 1
        if include:
            ptnkey = tuple(sorted(pattern))
            for read in pairs:
                read_patterns.append((read, ptnkey))
            nacounts += Counter(pattern)
            positions = set()
            for f, l in zip(firstpos, lastpos):
                if f < 0 or l < 0:
                    continue
                positions |= {
                    (pos, None, None)
                    for pos in range(f, l + 1)
                }
            nacounts += Counter(positions)

    # find > 10% mutations
    keepmuts = {}
    for (pos, na, refna), count in nacounts.items():
        if na is None:
            continue
        pcnt = count / nacounts[(pos, None, None)]
        if pcnt > MUTATION_PCNT_CUTOFF:
            keepmuts[(pos, na, refna)] = pcnt
    keepmutset = set(keepmuts.keys())

    for read, pattern in read_patterns:
        # remove <= 10% mutations
        pattern = tuple(sorted(set(pattern) & keepmutset))
        yield read, pattern


def load_patterns(pattern_csv, removes_csv):
    pattern_csv = csv.DictReader(pattern_csv)
    removes_csv = csv.reader(removes_csv)
    removes_indice = set()
    for _, indice in removes_csv:
        removes_indice |= {
            int(i) for i in indice.split(',') if i and i != 'None'
        }
    patterns = []
    for row in pattern_csv:
        index = int(row.pop('Index'))
        removed = index in removes_indice
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
            (index, samplename, ptid, pattern, count, pcnt, removed)
        )
    return patterns


def auto_is_removed_2(dist_interperson, dist_intraperson, dist_intrasample):
    outp_lt2 = np.count_nonzero(dist_interperson < 2)
    outp_lt5 = np.count_nonzero(dist_interperson < 5)
    if outp_lt2 == 0 and outp_lt5 < 2:
        return False  # should be kept
    if outp_lt2 > 3:
        return True  # should be removed
    if outp_lt2 == 0 and outp_lt5 > 50:
        pass
    if np.count_nonzero(dist_interperson < 2) > 0:
        return np.count_nonzero(dist_intraperson < 5) < 10
    else:
        return np.count_nonzero(dist_intraperson < 2) == 0


DIST_INTERPERSON_CUTOFF = 10
DIST_INTRAPERSON_CUTOFF = 14
DIST_INTRASAMPLE_CUTOFF = 12
MIN_FLOAT = sys.float_info.min


def calc_partial_score(dist_array, cutoff, sign=1):
    partials = []
    total = len(dist_array)
    maxdist = np.max(dist_array)
    for dist in range(0, max(cutoff + 1, maxdist)):
        # partial = (
        #     (np.e ** (sign * (cutoff - dist)) - 1) *
        #     np.count_nonzero(dist_array == dist) / total
        # )
        partial = (
            sign * (dist - cutoff) *
            np.count_nonzero(dist_array == dist) / total
        )
        if dist <= cutoff:
            partials.append(partial)
        else:
            partials[-1] += partial
        # partials.append((np.log10(
        #     np.count_nonzero(dist_array == dist) + MIN_FLOAT
        # ) - total_dist) * sign)
    # partials.append((
    #     np.log10(np.count_nonzero(dist_array > dist)
    #              + MIN_FLOAT) - total_dist
    # ) * sign * -1)
    # partials.append(
    #     -1 * sign * (np.count_nonzero(dist_array > dist) + MIN_FLOAT)
    # )
    return sum(partials), partials


def calc_score(dist_interperson, dist_intraperson, dist_intrasample):
    quantiles = np.array([0.01, 0.02, 0.05, 0.1])
    # inter-person
    outp = (np.quantile(dist_interperson, quantiles)
            if dist_interperson.size > 0 else np.array([-1] * quantiles.size))
    # intra-person
    inp = (np.quantile(dist_intraperson, quantiles)
           if dist_intraperson.size > 0 else np.array([-1] * quantiles.size))
    # intra-sample
    ins = (np.quantile(dist_intrasample, quantiles)
           if dist_intrasample.size > 0 else np.array([-1] * quantiles.size))
    return (np.max(
        np.concatenate((inp / outp, ins / outp))
    ), *outp, *inp, *ins)


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
        dist_interperson = np.array(sorted(dist_interperson))
        dist_intraperson = np.array(sorted(dist_intraperson))
        dist_intrasample = np.array(sorted(dist_intrasample))

        yield (
            calc_score(
                dist_interperson,
                dist_intraperson,
                dist_intrasample
            )
        )


@click.command()
@click.argument(
    'pattern_csv', type=click.File('r'))
@click.argument(
    'removes_csv', type=click.File('r'))
@click.argument(
    'output_csv', type=click.File('w'))
def main(pattern_csv, removes_csv, output_csv):
    patterns = load_patterns(pattern_csv, removes_csv)
    writer = csv.writer(output_csv)
    writer.writerow([
        'SampleName', 'Index', 'PtID', 'Pattern',
        'Count', 'Pcnt', 'IsRemoved',
        'Score',
        # *['OutP_eq{}'.format(d) for d in range(DIST_INTERPERSON_CUTOFF)],
        # 'OutP_gte{}'.format(DIST_INTERPERSON_CUTOFF),
        # *['InP_eq{}'.format(d) for d in range(DIST_INTRAPERSON_CUTOFF)],
        # 'InP_gte{}'.format(DIST_INTRAPERSON_CUTOFF),
        # *['InS_eq{}'.format(d) for d in range(DIST_INTRASAMPLE_CUTOFF)],
        # 'InS_gte{}'.format(DIST_INTRASAMPLE_CUTOFF),
        'OutP_5Pcnt',
        'InP_5Pcnt',
        'InS_5Pcnt',
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
    score_removed = []
    score_kept = []

    for row in iter_rows(patterns):
        (idx, samplename, ptid, pattern,
         count, pcnt, removed, *score) = row
        pattern = ', '.join(['{2}{0}{1}'.format(*m) for m in sorted(pattern)])
        if removed:
            score_removed.append(score[0])
        else:
            score_kept.append(score[0])
        writer.writerow([
            samplename, idx, ptid, pattern, count, pcnt, removed, *score
        ])
    click.echo('DIST_INTERPERSON_CUTOFF: {}'.format(DIST_INTERPERSON_CUTOFF))
    click.echo('DIST_INTRAPERSON_CUTOFF: {}'.format(DIST_INTRAPERSON_CUTOFF))
    click.echo('DIST_INTRASAMPLE_CUTOFF: {}'.format(DIST_INTRASAMPLE_CUTOFF))
    click.echo('Score Distribution:')
    quantiles = np.array([0, .1, .25, .5, .75, .9, 1])
    click.echo('  - Removed (0, 10, 25, 50, 75, 90, 100): ')
    click.echo('{} {} {} {} {} {} {}'.format(
        *np.quantile(score_removed, quantiles)))
    click.echo('  - Kept (0, 10, 25, 50, 75, 90, 100): ')
    click.echo('{} {} {} {} {} {} {}'.format(
        *np.quantile(score_kept, quantiles)))


if __name__ == '__main__':
    main()
