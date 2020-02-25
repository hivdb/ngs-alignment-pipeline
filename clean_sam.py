#! /usr/bin/env python

import os
import csv
import pysam
import click

import fastareader
from patternutils import iter_read_patterns

REF_CODON_OFFSET = 0

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 30))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 25))
MUTATION_PCNT_CUTOFF = 0.1

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


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
                _, lastrefprofile = fastareader.load(lastref)
            lastrefprofile = lastrefprofile['sequence']

            with pysam.AlignmentFile(sampath) as samfile:
                with pysam.AlignmentFile(cleaned_output, 'wb',
                                         template=samfile) as samout:
                    for read, pattern in iter_read_patterns(
                        sampath, initrefnas, lastrefprofile,
                        refbegin, refend, pos_offset,
                        MUTATION_PCNT_CUTOFF
                    ):
                        samplename = \
                            sampath.rsplit('/', 1)[-1].split('_', 1)[0]
                        if (samplename, pattern) not in pattern_keeps:
                            continue
                        samout.write(read)
                print(cleaned_output)


if __name__ == '__main__':
    main()
