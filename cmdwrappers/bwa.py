#! /usr/bin/env python

import os
from .base import docker_execute, refinit_func, align_func


@refinit_func('bwa')
def bwa_refinit(refseq):
    suffixes = ('amb', 'ann', 'bwt', 'pac', 'sa')
    if not all(os.path.isfile(
        os.path.extsep.join((refseq, suffix))
    ) for suffix in suffixes):
        # the reference index files are not previously generated
        refdir, refname = os.path.split(refseq)
        docker_execute([
            'bwa', 'index',
            '/shared/ref/{}'.format(refname)
        ], {
            refdir: {
                'bind': '/shared/ref',
                'mode': 'rw'
            }
        })


@align_func('bwa')
def bwa_align(refseq, fastq1, fastq2, sam, indir, outdir):
    refdir, refname = os.path.split(refseq)
    command = [
        'bwa', 'mem',
        '-o', '/shared/output/{}'.format(sam),
        '/shared/ref/{}'.format(refname),
        '/shared/input/{}'.format(fastq1),
    ]
    if fastq2:
        command.append('/shared/input/{}'.format(fastq2))
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
