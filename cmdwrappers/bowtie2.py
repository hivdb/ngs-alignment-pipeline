#! /usr/bin/env python

import os
import re
import click
from .base import docker_execute, refinit_func, align_func

BOWTIE2_ARGS = [
    '--local',
    '--threads', '3',
    '--rdg', '15,3',  # read gap open, extend penalties
    '--rfg', '15,3',  # reference gap open, extend penalties
    # '--ma', '2',     # match bonus
    # '--mp', '6',     # max mismatch penalty; lower qual = lower penalty
]


@refinit_func('bowtie2')
def bowtie2_refinit(refseq):
    refpath, ext = os.path.splitext(refseq)
    suffixes = ('1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2')
    if not all(os.path.isfile(
        os.path.extsep.join((refpath, suffix))
    ) for suffix in suffixes):
        # the reference index files are not previously generated
        refdir, refname = os.path.split(refpath)
        logs = docker_execute([
            'bowtie2-build',
            '/shared/ref/{}{}'.format(refname, ext),
            '/shared/ref/{}'.format(refname)
        ], {
            refdir: {
                'bind': '/shared/ref',
                'mode': 'rw'
            }
        })
        click.echo(logs.decode('U8'))


@align_func('bowtie2')
def bowtie2_align(refseq, fastq1, fastq2, sam, indir, outdir):
    refpath, _ = os.path.splitext(refseq)
    refdir, refname = os.path.split(refpath)
    command = [
        'bowtie2',
        *BOWTIE2_ARGS,
        '-x', '/shared/ref/{}'.format(refname),
        '-S', '/shared/output/{}'.format(sam),
        '-1', '/shared/input/{}'.format(fastq1),
    ]
    if fastq2:
        command.extend([
            '-2', '/shared/input/{}'.format(fastq2)
        ])
    logs = docker_execute(
        command,
        {
            refdir: {
                'bind': '/shared/ref',
                'mode': 'ro'
            },
            indir: {
                'bind': '/shared/input',
                'mode': 'ro'
            },
            outdir: {
                'bind': '/shared/output',
                'mode': 'rw'
            }
        }
    )
    overall_rate = re.search(
        r'\n(\d+\.\d+)% overall alignment rate',
        logs.decode('U8'))
    if overall_rate:
        overall_rate = float(overall_rate.group(1))
    else:
        overall_rate = -1
    with open(os.path.join(
        outdir, os.path.splitext(sam)[0] + '.log'
    ), 'wb') as fp:
        fp.write(logs)
    return {'overall_rate': overall_rate}
