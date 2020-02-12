#! /usr/bin/env python

import os
import re
from .base import docker_execute, refinit_func, align_func

SNAP_ARGS = [
    '-t', '3'
]


@refinit_func('snap')
def snap_refinit(refseq):
    refpath, ext = os.path.splitext(refseq)
    files = ('Genome', 'GenomeIndex', 'GenomeIndexHash', 'OverflowTable')
    if not all(os.path.isfile(
        os.path.join(refpath, fn)
    ) for fn in files):
        # the reference index files are not previously generated
        refdir, refname = os.path.split(refpath)
        docker_execute([
            'snap-aligner', 'index',
            '/shared/ref/{}{}'.format(refname, ext),
            '/shared/ref/{}'.format(refname)
        ], {
            refdir: {
                'bind': '/shared/ref',
                'mode': 'rw'
            }
        })


@align_func('snap')
def snap_align(refseq, fastq1, fastq2, sam, indir, outdir):
    refpath, _ = os.path.splitext(refseq)
    refdir, refname = os.path.split(refpath)
    subcommand = 'single'
    if fastq2:
        subcommand = 'paired'
    command = [
        'snap-aligner', subcommand,
        *SNAP_ARGS,
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
