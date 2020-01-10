import os
from .base import docker_execute


def stats(sam, outdir):
    logs = docker_execute([
        'samtools', 'stats',
        '/shared/output/{}'.format(sam)
    ], {
        outdir: {
            'bind': '/shared/output',
            'mode': 'rw'
        }
    })
    with open(os.path.join(
        outdir, os.path.splitext(sam)[0] + '.stats.txt'
    ), 'wb') as fp:
        fp.write(logs)
