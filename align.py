#! /usr/bin/env python

import os
import re
import click
import tempfile
import fastareader
from sam2fasta import sam2fasta
from itertools import combinations
from collections import defaultdict
from cmdwrappers import (
    samtools,
    get_programs, get_refinit, get_align
)

FILENAME_DELIMITERS = (' ', '_', '-')
PAIRED_FASTQ_MARKER = ('1', '2')


def find_paired_marker(text1, text2):
    diffcount = 0
    diffpos = -1
    for pos, (a, b) in enumerate(zip(text1, text2)):
        if diffcount > 1:
            return -1
        if a == b:
            continue
        if a not in PAIRED_FASTQ_MARKER or b not in PAIRED_FASTQ_MARKER:
            return -1
        diffcount += 1
        diffpos = pos
    return diffpos


def find_paired_fastq_patterns(filenames):
    """Smartly find paired FASTQ file patterns

    A valid filename pattern must meet:
    - use one of the valid delimiters (" ", "_" or "-") to separate the
      filename into different chunks
    - in one and only one chunk, a fixed position character changed from "1" to
      "2"

    Valid pair pattern examples:
      14258F_L001_R1_001.fastq.gz <-> 14258F_L001_R2_001.fastq.gz
      SampleExample_1.fastq <-> SampleExample_2.fastq

    Invalid pair pattern examples:
      SampleExample1.fastq <-> SampleExample2.fastq
      SampleExample_1.fastq <-> SampleExample_2.fastq.gz
      SampleExample_1.FASTQ.GZ <-> SampleExample_2.fastq.gz

    """
    patterns = defaultdict(list)
    for fn1, fn2 in combinations(filenames, 2):
        if len(fn1) != len(fn2):
            continue
        for delimiter in FILENAME_DELIMITERS:
            if delimiter not in fn1 or delimiter not in fn2:
                continue
            chunks1 = fn1.split(delimiter)
            chunks2 = fn2.split(delimiter)
            if len(chunks1) != len(chunks2):
                continue
            for reverse in range(2):
                diffcount = 0
                invalid = False
                diffoffset = -1
                if reverse:
                    chunks1.reverse()
                    chunks2.reverse()
                for n, (left, right) in enumerate(zip(chunks1, chunks2)):
                    if diffcount > 1:
                        invalid = True
                        break
                    if left == right:
                        continue
                    pos_paired_marker = find_paired_marker(left, right)
                    if pos_paired_marker < 0:
                        invalid = True
                        break
                    diffoffset = n
                    diffcount += 1
                if not invalid:
                    patterns[(
                        delimiter,
                        diffoffset,
                        pos_paired_marker,
                        reverse
                    )].append((fn1, fn2))
    covered = set()
    for pattern, pairs in sorted(
            patterns.items(), key=lambda p: (-len(p[1]), -p[0][3])):
        known = set()
        invalid = False
        for left, right in pairs:
            if left in covered or right in covered:
                # a pattern is invalid if the pairs is already matched
                # by a previous pattern
                invalid = True
                break

            if left in known or right in known:
                # a pattern is invalid if there's duplicate in pairs
                invalid = True
                break
            known.add(left)
            known.add(right)

        if not invalid:
            covered |= known
            yield pattern, pairs
    if len(filenames) > len(covered):
        remains = sorted(set(filenames) - covered)
        yield (None, -1, -1, -1), [(l, None) for l in remains]


def find_paired_fastqs(input_directory):
    for dirpath, _, filenames in os.walk(input_directory, followlinks=True):
        filenames = [
            fn for fn in filenames
            if fn[-6:].lower() == '.fastq'
            or fn[-9:].lower() == '.fastq.gz'
        ]
        yield from (
            (dirpath, fnpair, pattern)
            for pattern, fnpairs in find_paired_fastq_patterns(filenames)
            for fnpair in fnpairs
        )


def name_samfile(fnpair, pattern):
    delimiter, offset, _, reverse = pattern
    filename, _ = fnpair
    samfile = re.split(r'(?i)\.fastq', filename)[0]
    if reverse == -1:
        return samfile + '.sam'
    samfile = samfile.split(delimiter)
    if reverse:
        samfile.reverse()
    samfile = samfile[:offset] + samfile[offset + 1:]
    if reverse:
        samfile.reverse()
    return delimiter.join(samfile) + '.sam'


def replace_ext(filename, toext, fromext=None):
    if fromext:
        return filename[-len(fromext):] + toext
    else:
        return os.path.splitext(filename)[0] + toext


def fasta_first_sequence(filename):
    with open(filename) as fp:
        return fastareader.load(fp)[0]['sequence']


def calc_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise RuntimeError('sequences are not the same size')
    return sum(na1 != na2 for na1, na2 in zip(seq1, seq2)) / len(seq1)


@click.command()
@click.argument(
    'input_directory',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.argument(
    'output_directory',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.option(
    '-p', '--program',
    required=True,
    type=click.Choice(get_programs()))
@click.option(
    '-r', '--reference',
    required=True,
    type=click.Path(exists=True, file_okay=True,
                    dir_okay=False, resolve_path=True))
@click.option('-i', '--iterative', is_flag=True)
@click.option('--max-iterative-circles', default=20, type=int)
def main(
    input_directory, output_directory,
    program, reference, iterative,
    max_iterative_circles
):
    refinit = get_refinit(program)
    align = get_align(program)
    if not iterative:
        max_iterative_circles = 1
    for dirpath, fnpair, pattern in find_paired_fastqs(input_directory):
        samfile = name_samfile(fnpair, pattern)
        refseq = reference
        prevrefnas = fasta_first_sequence(refseq)
        all_prevrefnas = {prevrefnas}
        with tempfile.TemporaryDirectory(dir="/tmp") as tmpdir:
            for i in range(max_iterative_circles):
                suffix = '.{}'.format(
                    i if i + 1 < max_iterative_circles else 'final')
                refinit(refseq)
                align(
                    refseq, *fnpair, samfile,
                    dirpath, output_directory)
                samtools.stats(samfile, output_directory)
                # bamfile = replace_ext(samfile, suffix + '.bam')
                fastafile = replace_ext(samfile, suffix + '.fas')
                refnas = sam2fasta(
                    prevrefnas, os.path.join(output_directory, samfile))
                if refnas in all_prevrefnas:
                    refseq = os.path.join(output_directory, fastafile)
                else:
                    refseq = os.path.join(tmpdir, fastafile)
                with open(refseq, 'w') as fp:
                    fp.write('>Consensus\n')
                    fp.write(refnas)
                if refnas in all_prevrefnas:
                    click.echo('{} round {}: completed'.format(samfile, i + 1))
                    break
                click.echo(
                    '{} round {}: consensus distance={:.4f}%'
                    .format(samfile, i + 1,
                            calc_distance(prevrefnas, refnas) * 100))
                all_prevrefnas.add(refnas)
                prevrefnas = refnas


if __name__ == '__main__':
    main()
