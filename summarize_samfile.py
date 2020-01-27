#! /usr/bin/env python

import os
import csv
import click

CSV_HEADER = [
    'Sample name',
    'Overall alignment rate (%)',
    'Sequences',
    'Reads mapped',
    'Reads unmapped',
    '',
    'Raw total sequences',
    'Filtered sequences',
    'Is sorted',
    '1st fragments',
    'Last fragments',
    'Reads mapped and paired',
    'Reads properly paired',
    'Reads paired',
    'Reads duplicated',
    'Reads MQ0',
    'Reads QC failed',
    'Non-primary alignments',
    'Total length',
    'Total first fragment length',
    'Total last fragment length',
    'Bases mapped',
    'Bases mapped (cigar)',
    'Bases trimmed',
    'Bases duplicated',
    'Mismatches',
    'Error rate',
    'Average length',
    'Average first fragment length',
    'Average last fragment length',
    'Maximum length',
    'Maximum first fragment length',
    'Maximum last fragment length',
    'Average quality',
    'Insert size average',
    'Insert size standard deviation',
    'Inward oriented pairs',
    'Outward oriented pairs',
    'Pairs with other orientation',
    'Pairs on different chromosomes',
    'Percentage of properly paired reads (%)',
]


@click.command()
@click.argument(
    'input_directory',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.argument('output_csv', type=click.File('w'))
def main(input_directory, output_csv):
    writer = csv.DictWriter(output_csv, CSV_HEADER)
    writer.writeheader()
    for indir, _, filenames in os.walk(input_directory):
        for fn in filenames:
            if not fn.endswith('.stats.txt'):
                continue
            with open(os.path.join(indir, fn)) as fp:
                row = {
                    'Sample name': fn[:-10],
                    '': ''
                }
                for line in fp:
                    if line.startswith('SN'):
                        _, key, value, *_ = line.split('\t')
                        key = key[0].upper() + key[1:-1]
                        row[key] = value.strip()
                    elif line.startswith('FFQ'):
                        break
                row['Overall alignment rate (%)'] = '{:.2f}'.format(
                    100 * float(row['Reads mapped']) / float(row['Sequences'])
                )
                writer.writerow(row)


if __name__ == '__main__':
    main()
