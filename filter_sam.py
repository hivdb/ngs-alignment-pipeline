#! /usr/bin/env python

import os
import csv
import pysam
import click
from statistics import mean
from collections import defaultdict, Counter, OrderedDict

import fastareader

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
    if mean(qua) < OVERALL_QUALITY_CUTOFF:
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


def load_pattern_keeps(pattern_csv, removes_csv):
    pattern_csv = csv.DictReader(pattern_csv)
    removes_csv = csv.reader(removes_csv)
    removes_indice = set()
    for _, indice in removes_csv:
        removes_indice |= {
            int(i) for i in indice.split(',') if i and i != 'None'
        }
    keeps = set()
    for row in pattern_csv:
        index = int(row.pop('Index'))
        if index in removes_indice:
            continue
        samplename = row.pop('SampleName')
        row.pop('Count')
        row.pop('Pcnt')
        pattern = set()
        for pos, na in row.items():
            if na in ('-', '.'):
                continue
            ref = pos[0]
            pos = int(pos[1:])
            pattern.add((pos, na, ref))
        keeps.add((samplename, tuple(sorted(pattern))))
    return keeps


@click.command()
@click.argument(
    'initref', type=click.File('r'))
@click.argument(
    'samfolder', type=click.Path(dir_okay=True, exists=True))
@click.argument(
    'pattern_csv', type=click.File('r'))
@click.argument(
    'removes_csv', type=click.File('r'))
@click.option('--ref-range', default='0-99999')
@click.option('--pos-offset', type=int, default=0)
def main(samfolder, initref, pattern_csv, removes_csv, ref_range, pos_offset):
    refbegin, refend = [int(r) for r in ref_range.split('-', 1)]
    initrefnas, = fastareader.load(initref)
    initrefnas = initrefnas['sequence']
    pattern_keeps = load_pattern_keeps(pattern_csv, removes_csv)
    for basedir, _, filenames in os.walk(samfolder):
        for sampath in filenames:
            if sampath[-4:].lower() not in ('.sam', '.bam'):
                continue
            if sampath.endswith('.cleaned.bam'):
                continue
            sampath = os.path.join(basedir, sampath)
            lastref = os.path.splitext(sampath)[0] + '.lastref.fas'
            # mutcoroutput = os.path.splitext(sampath)[0] + '.mutcor.csv'
            cleaned_output = os.path.splitext(sampath)[0] + '.cleaned.bam'
            with open(lastref) as lastref:
                _, alnprofile = fastareader.load(lastref)
            # refnas = refnas['sequence']
            alnprofile = alnprofile['sequence']
            alnprofile = attach_initref_pos(alnprofile)
            with pysam.AlignmentFile(sampath) as samfile:
                with pysam.AlignmentFile(cleaned_output, 'wb',
                                         template=samfile) as samout:
                    for read, pattern in find_patterns(
                        sampath, initrefnas, alnprofile,
                        refbegin, refend, pos_offset
                    ):
                        samplename = \
                            sampath.rsplit('/', 1)[-1].split('_', 1)[0]
                        if (samplename, pattern) not in pattern_keeps:
                            continue
                        samout.write(read)
                print(cleaned_output)


if __name__ == '__main__':
    main()
