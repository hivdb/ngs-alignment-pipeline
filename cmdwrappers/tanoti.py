#! /usr/bin/env python

import os
from .base import docker_execute, refinit_func, align_func

TANOTI_ARGS = [
    '-m', '50',  # minimal match % of a read
    '-u', '1',   # include unmapped reads in the output
    '-t', '1',   # keep temporary files
    '-P', '8',   # number of parallel BLAST search
]


@refinit_func('tanoti')
def tanoti_refinit(refseq):
    pass


@align_func('tanoti')
def tanoti_align(refseq, fastq1, fastq2, sam, indir, outdir):
    refdir, refname = os.path.split(refseq)
    command = [
        'tanoti-wrapper',
        *TANOTI_ARGS,
        '-r', '/shared/ref/{}'.format(refname),
        '-i', '/shared/input/{}'.format(fastq1),
    ]
    if fastq2:
        command.extend([
            '/shared/input/{}'.format(fastq2),
            '-p', '1'
        ])
    command.extend([
        '-o', '/shared/output/{}'.format(sam)
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
    overall_rate = -1.
    with open(os.path.join(
        outdir, os.path.splitext(sam)[0] + '.log'
    ), 'wb') as fp:
        fp.write(logs)
    return {'overall_rate': overall_rate}
