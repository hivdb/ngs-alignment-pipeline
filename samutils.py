#! /usr/bin/env python

import os
import pysam
from statistics import mean
from collections import defaultdict, OrderedDict, Counter

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 30))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 25))

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


def extract_nas(read):
    seq, qua, aligned_pairs = (
        read.query_sequence,
        read.query_qualities,
        read.get_aligned_pairs(False)
    )

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
        result_nas.append((pos, ''.join(n for n, _ in nqs), meanq))
    if result_nas:
        lastpos, lastna, q = result_nas[-1]
        # remove insertion at the end of sequence read
        result_nas[-1] = (lastpos, lastna[0], q)
    return result_nas, err


def iter_paired_reads(samfile):
    paired_reads = OrderedDict()
    with pysam.AlignmentFile(samfile, 'rb') as samfile:
        for idx, read in enumerate(samfile.fetch()):
            name = read.query_name
            paired_reads.setdefault(name, []).append(read)
    for header, pair in paired_reads.items():
        if len(pair) > 2:
            raise RuntimeError(
                'Malformed SAM file: too many reads in one pair')
        yield header, pair


def iter_posnas(all_paired_reads,
                site_quality_cutoff=SITE_QUALITY_CUTOFF):

    for header, pair in all_paired_reads:
        posnas = defaultdict(Counter)
        for read in pair:
            results, err = extract_nas(read)
            if err & ERR_TOO_SHORT or err & ERR_LOW_QUAL:
                continue
            for refpos, nas, q in results:
                if q < site_quality_cutoff:
                    continue
                posnas[refpos][nas] = q
        if posnas:
            yield header, [
                # pos, nas, qual
                (pos, *counter.most_common(1)[0])
                for pos, counter in sorted(posnas.items())
            ]
