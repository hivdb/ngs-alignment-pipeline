#! /usr/bin/env python

from collections import defaultdict, Counter

from samreader import iter_paired_reads


def sam2fasta(refnas, refalnprofile, sampath):
    nafreqs = defaultdict(Counter)
    for _, allnas in iter_paired_reads(sampath):
        for refpos, nas, _ in allnas:
            for i, na in enumerate(nas):
                nafreqs[(refpos, i)][na] += 1

    resultseq = []
    alnprofile = []
    for refpos0, (refna, refp) in enumerate(zip(refnas, refalnprofile)):
        refpos = refpos0 + 1
        if (refpos, 0) in nafreqs:
            counter = nafreqs[(refpos, 0)]
            total = sum(counter.values())
            (na, _), = counter.most_common(1)
            resultseq.append(refna if na == '-' else na)
            alnprofile.append(
                '+' if refp == '+' else '-' if na == '-' else ':')
            i = 1
            while (refpos, i) in nafreqs:
                counter = nafreqs[(refpos, i)]
                i += 1
                (na, count), = counter.most_common(1)
                if count * 2 >= total:
                    resultseq.append(na)
                    alnprofile.append('+')
                else:
                    break
        else:
            resultseq.append(refna)
            alnprofile.append('.')
    return ''.join(resultseq), ''.join(alnprofile)
