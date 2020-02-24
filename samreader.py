#! /usr/bin/env python

import os
import pysam
from statistics import mean
from multiprocessing import Process, Queue
from collections import defaultdict, Counter

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 30))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 25))

NUM_PROCESSES = int(os.environ.get('NTHREADS', 3))
INPUT_QUEUE = Queue(NUM_PROCESSES)
OUTPUT_QUEUE = Queue()

PRODUCER_CAPACITY = 50000

CHUNKSIZE = 500

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


def get_na_counts(seq, qua, aligned_pairs, idx, header):

    # pre-filter
    err = ERR_OK
    if len(seq) < LENGTH_CUTOFF:
        err |= ERR_TOO_SHORT
    if mean(qua) < OVERALL_QUALITY_CUTOFF:
        err |= ERR_LOW_QUAL
    if err:
        return (header, None), err

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
        meanq = mean(qs)
        result_nas.append((pos, ''.join(n for n, _ in nqs), meanq))
    if result_nas:
        lastpos, lastna, q = result_nas[-1]
        # remove insertion at the end of sequence read
        result_nas[-1] = (lastpos, lastna[0], q)
    return (header, result_nas), err


def reads_consumer():
    while True:
        chunk = INPUT_QUEUE.get()
        out_chunk = []
        for args in chunk:
            out_chunk.append(get_na_counts(*args))
        OUTPUT_QUEUE.put(out_chunk)


def reads_producer(filename, offset, producer_capacity):
    with pysam.AlignmentFile(filename, 'rb') as samfile:
        chunk = []
        limit = offset + producer_capacity
        for idx, read in enumerate(samfile.fetch()):
            if idx < offset:
                continue
            if idx >= limit:
                break
            seq, qua, aligned_pairs, name = (read.query_sequence,
                                             read.query_qualities,
                                             read.get_aligned_pairs(False),
                                             read.query_name)
            if len(chunk) == CHUNKSIZE:
                INPUT_QUEUE.put(chunk)
                chunk = []
            chunk.append((seq, qua, aligned_pairs, idx + offset, name))
        if chunk:
            INPUT_QUEUE.put(chunk)


def create_producers(sampath, totalreads, producer_capacity):
    offset = 0
    processes = []
    while offset < totalreads:
        producer = Process(target=reads_producer, args=(
            sampath, offset, producer_capacity))
        producer.daemon = True
        processes.append(producer)
        producer.start()
        offset += producer_capacity
    return processes


def create_consumers(num_processes):
    consumers = []
    for _ in range(num_processes):
        consumer = Process(target=reads_consumer)
        consumer.daemon = True
        consumers.append(consumer)
        consumer.start()
    return consumers


def iter_paired_reads(
    sampath,
    site_quality_cutoff=SITE_QUALITY_CUTOFF,
    producer_capacity=PRODUCER_CAPACITY
):
    with pysam.AlignmentFile(sampath, 'rb') as samfile:
        # set until_eof=False to exclude unmapped reads
        totalreads = samfile.count(until_eof=False)

    producers = create_producers(sampath, totalreads, PRODUCER_CAPACITY)

    num_finished = 0
    paired_results = defaultdict(list)
    while num_finished < totalreads:
        out_chunk = OUTPUT_QUEUE.get()
        for (header, my_results), my_err in out_chunk:
            num_finished += 1
            paired_results[header].append((my_results, my_err))
    [p.terminate() for p in producers]

    for header, pairs in paired_results.items():
        if len(pairs) > 2:
            raise RuntimeError(
                'Malformed SAM file: too many reads in one pair')
        allnas = defaultdict(Counter)
        for j, (results, err) in enumerate(pairs):
            if err & ERR_TOO_SHORT or err & ERR_LOW_QUAL:
                continue
            for refpos, nas, q in results:
                if q < site_quality_cutoff:
                    continue
                allnas[refpos][nas] = q
        if allnas:
            yield header, [
                # pos, nas, qual
                (pos, *counter.most_common(1)[0])
                for pos, counter in sorted(allnas.items())
            ]


consumers = create_consumers(NUM_PROCESSES)
