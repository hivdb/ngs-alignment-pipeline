import re
import ssw
from collections import defaultdict, Counter
from samutils import iter_paired_reads, iter_posnas


CIGAR_PATTERN = re.compile(r'(\d+)([MNDI])')


def str_index_all(text, char):
    start = 0
    indices = []
    while True:
        try:
            idx = text.index(char, start)
            indices.append(idx)
            start = idx + 1
        except ValueError:
            break
    return indices


def realign(initref, lastref):
    """
    Realign lastref to initref

    The purpose is to generate a better alignment for restoring
    initref position numbers in multi-alignment results.
    """
    cons_delpos = set(str_index_all(lastref, '-'))
    lastref = lastref.replace('-', '')
    aligner = ssw.Aligner()
    aln = aligner.align(reference=initref, query=lastref)
    if aln.query_begin > 0 or aln.reference_begin > 0:
        # Why?
        raise RuntimeError('Consensus misaligned')
    initref, _, lastref = aln.alignment
    ref_pos0 = -1
    alnprofile = []
    for refna, consna in zip(initref, lastref):
        if refna != '-':
            ref_pos0 += 1

        # The re-alignment could place deletions in slightly different
        # place comparing to the original alignment. Four cases need
        # to be considered:
        #   - +n+o: deletion presents in both alignments:
        #           no extra handling needed
        #   - +n-o: deletion only presents in new alignment:
        #           new deletion should be added to the multi-alignment
        #   - -n+o: deletion only presents in old alignment
        #           multi-alignment NAs mapped to the old deletion
        #           should be removed to previous refpos
        #   - -n-o: no deletion:
        #           no extra handling needed
        newdel = consna == '-'
        olddel = len(alnprofile) in cons_delpos
        while not newdel and olddel:  # -n+o
            # mismatched old deletions
            prev_refpos0, _ = alnprofile[-1]
            alnprofile.append((prev_refpos0, -1))
            olddel = len(alnprofile) in cons_delpos
        else:
            if newdel:  # +n+o/+n-o
                # new deletions
                prev_refpos0, count = alnprofile[-1]
                alnprofile[-1] = (prev_refpos0, count + 1)
            else:  # -n-o
                # agree, non deletion
                alnprofile.append((ref_pos0, 0))
    return alnprofile


def map_posnas_to_initref(posnas, initrefnas, alnprofile):
    """
    Re-map an iteratively aligned reads to the original reference

    The position numbers of original reference (initref) was lost during
    the iterative alignment. This generator function restores the position
    numbers by using information from alignment profile list (alnprofile).

    This generator produces a tuple of three elements:
      - `position`: The original reference position minus `pos_offset`.
      - `nas`: Nucleic acid notations of this position. The insertion gap is
        represented in higher numbers (>1) of `nas`.
      - `refna`: Nucleic acid notation of the orig. reference at this positon.
    """
    initrefpos_map = defaultdict(list)
    for refpos, nas, _ in posnas:
        refpos0 = refpos - 1
        initrefpos0, modifier = alnprofile[refpos0]
        initrefpos_map[initrefpos0].append(nas)
        while modifier > 0:
            initrefpos0 += 1
            modifier -= 1
            initrefpos_map[initrefpos0].append('-')
    for initrefpos0, nas in sorted(initrefpos_map.items()):
        nas = ''.join(nas)
        if len(nas) > 1:
            # remove next-to-insertion deletions
            nas = nas.replace('-', '')
        initrefpos = initrefpos0 + 1
        refna = initrefnas[initrefpos0]
        yield initrefpos, nas, refna


def trim_and_offset_posnas(posnas, refbegin, refend, ref_offset):
    """
    Trim input `posnas` and offset the position by `ref_offset
    """
    for refpos, nas, refna in posnas:
        if refpos < refbegin or refpos > refend:
            continue
        rel_refpos = refpos - ref_offset
        yield rel_refpos, nas, refna


def replace_indel_notations(posnas):
    """
    Replace the indel notations from long NAs and "-" to "ins" and "del"
    """
    for refpos, nas, refna in posnas:
        if len(nas) > 1:
            yield refpos, 'ins', refna
        else:
            yield refpos, nas.replace('-', 'del'), refna


def extract_mutations(posnas):
    """
    Extract mutations (diff from refna)

    Here we assume that refna doesn't have indels
    """
    return ((p, n, r) for p, n, r in posnas if n != r)


def extract_positions(posnas):
    """Extract positions"""
    return (p for p, _, _ in posnas)


def adjust_posnas(all_posnas, initrefnas, alnprofile,
                  refbegin, refend, pos_offset):
    """Apply pre-process steps for header_posnas"""
    for header, posnas in all_posnas:

        # pre-process posnas
        posnas = map_posnas_to_initref(posnas, initrefnas, alnprofile)
        posnas = trim_and_offset_posnas(posnas, refbegin, refend, pos_offset)
        posnas = set(replace_indel_notations(posnas))

        yield header, posnas


def iter_patterns(all_posnas):
    for header, posnas in all_posnas:

        # extract mutations & positions
        mutations = set(extract_mutations(posnas))
        coverage = set(extract_positions(posnas))

        yield (tuple(sorted(mutations)),
               tuple(sorted(coverage)),
               bool(posnas))


def stat_posnas(all_posnas):
    nacounts = Counter()
    for (_, posnas) in all_posnas:
        nacounts += Counter(posnas)
        nacounts += Counter({(pos, None, refna)
                             for pos, _, refna in posnas})
    return nacounts


def stat_patterns(patterns):
    pattern_counts = Counter()
    for pattern, cov, posna_not_empty in patterns:
        if posna_not_empty:
            pattern_counts[(pattern, cov)] += 1
    return pattern_counts


def find_usual_muts(patterns, nacounts, mutation_pcnt_cutoff):
    mutcounts = Counter()
    for pattern, cov, posna_not_empty in patterns:
        if posna_not_empty:
            # mutation counter ++
            mutcounts += Counter(pattern)
    # find out mutations about mutation_pcnt_cutoff
    keepmuts = {}
    for (pos, na, refna), count in mutcounts.items():
        pcnt = count / nacounts[(pos, None, refna)]
        if pcnt >= mutation_pcnt_cutoff:
            keepmuts[(pos, na, refna)] = pcnt
    return keepmuts


def patterns_without_unusual_muts(pattern_counts, usualmuts):
    for (old_pattern, coverage), count in pattern_counts.items():
        new_pattern = tuple(sorted(set(old_pattern) & usualmuts))
        yield new_pattern, coverage, count


def _internal_find_patterns(sampath, initrefnas, alnprofile,
                            refbegin, refend, pos_offset,
                            mutation_pcnt_cutoff):
    all_paired_reads = list(iter_paired_reads(sampath))
    all_posnas = iter_posnas(all_paired_reads)
    all_posnas = list(
        adjust_posnas(all_posnas, initrefnas, alnprofile,
                      refbegin, refend, pos_offset)
    )
    patterns = list(iter_patterns(all_posnas))

    nacounts = stat_posnas(all_posnas)
    pattern_counts = stat_patterns(patterns)

    # remove mutations below MUTATION_PCNT_CUTOFF from patterns and re-count
    keepmuts = find_usual_muts(patterns, nacounts, mutation_pcnt_cutoff)
    keepmutset = set(keepmuts.keys())
    final_patterns = patterns_without_unusual_muts(pattern_counts, keepmutset)

    return all_paired_reads, final_patterns, keepmuts, nacounts


def iter_read_patterns(sampath, initrefnas, alnprofile,
                       refbegin, refend, pos_offset,
                       mutation_pcnt_cutoff):
    all_paired_reads, patterns, _, _ = _internal_find_patterns(
        sampath, initrefnas, alnprofile,
        refbegin, refend, pos_offset,
        mutation_pcnt_cutoff)

    for (_, pair), (pattern, _, _) in zip(all_paired_reads, patterns):
        for read in pair:
            yield read, pattern


def find_patterns(sampath, initrefnas, alnprofile,
                  refbegin, refend, pos_offset,
                  mutation_pcnt_cutoff):
    _, patterns, keepmuts, nacounts = _internal_find_patterns(
        sampath, initrefnas, alnprofile,
        refbegin, refend, pos_offset,
        mutation_pcnt_cutoff)

    pattern_counts = Counter()
    for pattern, cov, count in patterns:
        pattern_counts[(pattern, cov)] += count
    return pattern_counts, keepmuts, nacounts


def filter_patterns(patterns, muts,
                    pattern_pcnt_idv_cutoff,
                    pattern_pcnt_acc_cutoff,
                    missing_position_threshold,
                    remove_partials=True):
    """Apply filter rules to patterns generated by `find_patterns`

    Rules:

      - MISSING_POSITION_THRESHOLD: Without remove_partials=True, the
        function allows to keep a partial pattern, if its position coverage
        is above this threshold. The pattern will be furtherly judged by
        the two below rules.
      - PATTERN_PCNT_ACC_CUTOFF: Sorted patterns by percent decendingly,
        then retain the top PATTERN_PCNT_ACC_CUTOFF (0-1.0) of patterns.
      - PATTERN_PCNT_IDV_CUTOFF: Patterns above this threshold are all
        kept no matter to PATTERN_PCNT_ACC_CUTOFF set.

    """
    collapsed = Counter()
    total = sum(patterns.values())
    idv_cutoff = total * pattern_pcnt_idv_cutoff
    acc_cutoff = total * pattern_pcnt_acc_cutoff
    mutpositions = {p for p, _, _ in muts}
    for (pattern, coverage), count in patterns.items():
        coverage = set(coverage)
        if remove_partials:
            # Partial covered pattern is not allowed
            # if specified remove_partials=True.
            if not coverage.issupperset(mutpositions):
                continue
        else:
            # Partial covered pattern is allowed.
            # Populate missed positions into the pattern
            pattern = set(pattern)
            missed = set()
            for pos, _, ref in muts:
                if pos not in coverage:
                    pattern.add((pos, '.', ref))
                    missed.add(pos)

            # check and skip if the number of missed
            # positions exceeds threshold
            if len(missed) > missing_position_threshold:
                continue
            pattern = tuple(sorted(pattern))

        collapsed[pattern] += count

    results = []
    countsum = 0
    for pattern, count in collapsed.most_common():
        if countsum > acc_cutoff and count < idv_cutoff:
            break
        results.append([pattern, count, count / total])
        countsum += count
    return results
