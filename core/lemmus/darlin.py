"""Lemmus：DARLIN 流式提取、共识与 clone 交付。"""
from __future__ import annotations
from .io import project_temporary_parent, atomic_output

import argparse
import csv
import gzip
import hashlib
import io
import os
import re
import sys
import time
import tempfile
from contextlib import ExitStack
from itertools import islice
from pathlib import Path
import numpy as np
import pandas as pd
from .config import (
    DEFAULT_EXTRACTION,
    DEFAULT_QC,
    select_config,
)
from .inputs import (_barcode_map_digest)


CA_5prime = 'AGCTGTACAAGTAAGCGGC'


CA_5prime_full = 'GAGCTGTACAAGTAAGCGGC'


CA_3prime = 'AGAATTCTAACTAGAGCTCGCTGATCAGCCTCGACTGTGCCTTCT'


CA_CARLIN = 'CGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCGACTATGGAGTCGCGACTGTACGCACACGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGAGTCGAGACGCTGACGATATGGAGTCGATACGTAGCACGCAGACGATGGGAGCT'


TA_5prime = 'TCGGTACCTCGCGAA'


TA_5prime_full = 'GCTCGGTACCTCGCGAA'


TA_3prime = 'GTCTTGTCGGTGCCTTCTAGTT'


TA_CARLIN = 'TCGCCGGAGTCGAGACGCTGACGATATGGAGTCGACACGACTCGCGCATACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGATACGATACGCGCACGCTATGGAGTCGCGACTGTACGCACACGCGATGGAGTCGACTGCACGACAGTCGACGATGGAGTCGATACGTAGCACGCAGACGATGGGAGCGAGAGCGCGCTCGTCGACTATGGA'


RA_5prime = 'GTACAAGTAAAGCGGCC'


RA_5prime_full = 'ATGTACAAGTAAAGCGGCC'


RA_3prime = 'GTCTGCTGTGTGCCTTCTAGTT'


RA_CARLIN = 'GCGCCGGCGAGCGCTATGAGCGACTATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGATACGATACGCGCACGCTATGGAGTCGACTGCACGACAGTCGACGATGGAGTCGATACGTAGCACGCAGACGATGGGAGCGAGTCGAGACGCTGACGATATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGCGACTGTACGCACACGCGATGGAGTCGAGAGCGCGCTCGTCGACTATGGA'


def _template_parts(template, use_short_3prime=False):
    """Return primer and reference sequence parts for one DARLIN template."""
    if template.startswith('cCARLIN'):
        appendix = 'TGGGAGCT' if use_short_3prime else ''
        seq_3prime = 'TGGGAGCTAGAA' if use_short_3prime else CA_3prime
        return (CA_5prime, CA_5prime_full, seq_3prime, appendix, CA_CARLIN)
    if template.startswith('Tigre'):
        appendix = 'TGGA' if use_short_3prime else ''
        seq_3prime = 'TGGAGTCTTG' if use_short_3prime else TA_3prime
        return (TA_5prime, TA_5prime_full, seq_3prime, appendix, TA_CARLIN)
    if template.startswith('Rosa'):
        appendix = 'TGGA' if use_short_3prime else ''
        seq_3prime = 'TGGAGTCTGC' if use_short_3prime else RA_3prime
        return (RA_5prime, RA_5prime_full, seq_3prime, appendix, RA_CARLIN)
    raise ValueError("template must start with {'cCARLIN','Tigre','Rosa'}")


def check_editing(df, template, progress=None):
    """Mark clone IDs that differ from the unedited template reference."""
    if progress is None:
        progress = print
    CARLIN_seq = _template_parts(template)[4]
    required_columns = {'cell_id', 'clone_id'}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(f'DARLIN editing input lacks required columns: {sorted(missing_columns)}')
    if df.empty:
        df['edited'] = pd.Series(False, index=df.index, dtype=bool)
        progress('Editing efficiency: nan (no valid DARLIN reads)')
        return df
    df['edited'] = ~df['clone_id'].astype(str).str.startswith(CARLIN_seq)
    editing_efficiency = df.groupby('cell_id')['edited'].mean().mean()
    progress(f'Editing efficiency: {editing_efficiency:.3f}')
    return df


def find_approximate_subsequence(sequence, target, max_mismatches=0):
    """Return the best target match index when mismatches are within the limit.

    Seed-bounded substitution search: a window carrying at most
    ``max_mismatches`` substitutions must contain at least one error-free
    chunk when the target is partitioned into ``max_mismatches + 1`` disjoint
    chunks, so exact chunk hits bound the candidate starts instead of paying
    an O(read_length x target_length) scan per read. Ties keep the earliest
    index, matching the historical exhaustive scan.
    """
    max_mismatches = int(max_mismatches)
    if max_mismatches <= 0 or not sequence or not target or len(sequence) < len(target):
        return (None, None)
    n = len(target)

    def best_window(candidate_starts):
        best_index = None
        best_mismatches = n + 1
        for start in sorted(candidate_starts):
            window = sequence[start:start + n]
            mismatches = sum(a != b for a, b in zip(window, target, strict=True))
            if mismatches < best_mismatches:
                best_index, best_mismatches = start, mismatches
                if mismatches == 0:
                    break
        return (best_index, best_mismatches) if best_mismatches <= max_mismatches else (None, None)

    chunk_count = max_mismatches + 1
    chunk_length = n // chunk_count
    if chunk_length < 1:
        # Target too short to seed; the exhaustive scan is already cheap here.
        return best_window(range(len(sequence) - n + 1))
    remainder = n - chunk_length * chunk_count
    candidate_starts = set()
    offset = 0
    for chunk_index in range(chunk_count):
        length = chunk_length + (1 if chunk_index < remainder else 0)
        chunk = target[offset:offset + length]
        pos = sequence.find(chunk)
        while pos >= 0:
            start = pos - offset
            if 0 <= start <= len(sequence) - n:
                candidate_starts.add(start)
            pos = sequence.find(chunk, pos + 1)
        offset += length
    return best_window(candidate_starts)


def _bounded_edit_distance(left, right, max_distance=1):
    """Return Levenshtein distance up to ``max_distance``; otherwise max+1."""
    left = str(left)
    right = str(right)
    max_distance = int(max_distance)
    if abs(len(left) - len(right)) > max_distance:
        return max_distance + 1
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        current = [i]
        row_min = current[0]
        for j, b in enumerate(right, 1):
            value = min(
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + (a != b),
            )
            current.append(value)
            row_min = min(row_min, value)
        if row_min > max_distance:
            return max_distance + 1
        previous = current
    return previous[-1]


def find_indel_aware_subsequence(sequence, target, max_edits=1, downstream_reference=None):
    """Locate a conserved 5' anchor with a seed-bounded edit-distance-1 search.

    For edit distance <=1, at least one of two non-overlapping target halves must
    match exactly.  We use those exact seeds to infer a very small set of candidate
    starts, then validate only target-width / +/-1 windows.  This avoids the
    previous O(read_length * primer_length) Python scan on every exact-match miss.
    Boundary indels within the final three primer bases remain fail-closed so a
    real DARLIN allele base is never consumed as primer sequence.
    """
    sequence = str(sequence).upper()
    target = str(target).upper()
    exact = sequence.find(target)
    if exact >= 0:
        return (exact, exact + len(target), 0, 'exact')
    max_edits = int(max_edits)
    if max_edits < 1 or not sequence or not target:
        return (None, None, None, 'none')

    n = len(target)
    protected_tail = 3
    split = max(1, n // 2)
    seeds = ((target[:split], 0), (target[split:], split))
    candidate_starts = set()
    for seed, offset in seeds:
        pos = sequence.find(seed)
        while pos >= 0:
            # A one-base insertion/deletion before the seed shifts the inferred
            # primer start by at most one base.
            inferred = pos - offset
            for delta in (-1, 0, 1):
                start = inferred + delta
                if 0 <= start < len(sequence):
                    candidate_starts.add(start)
            pos = sequence.find(seed, pos + 1)

    # If neither half is present exactly, edit distance 1 is impossible.
    if not candidate_starts:
        return (None, None, None, 'none')

    substitution_hits = []
    deletion_hits = []   # read is one base shorter than primer
    insertion_hits = []  # read is one base longer than primer
    for start in sorted(candidate_starts):
        # substitution-only: same width, exactly one mismatch
        if start + n <= len(sequence):
            window = sequence[start:start+n]
            mismatches = sum(a != b for a, b in zip(window, target, strict=True))
            if mismatches == 1:
                substitution_hits.append((start, start+n, 1, 'substitution'))

        # deletion in read relative to primer.  Validate the deleted target
        # position and reject boundary deletions near the editable array.
        if n > 1 and start + n - 1 <= len(sequence):
            window = sequence[start:start+n-1]
            for deleted_pos in range(0, max(0, n - protected_tail)):
                if target[:deleted_pos] + target[deleted_pos+1:] == window:
                    deletion_hits.append((start, start+n-1, 1, 'deletion'))
                    break

        # insertion in read relative to primer.  The inserted base must be away
        # from the protected primer/allele boundary.
        if start + n + 1 <= len(sequence):
            window = sequence[start:start+n+1]
            for inserted_pos in range(0, max(0, n - protected_tail + 1)):
                if window[:inserted_pos] + window[inserted_pos+1:] == target:
                    insertion_hits.append((start, start+n+1, 1, 'insertion'))
                    break

    def unique_hit(hits):
        unique = {(x[0], x[1], x[3]): x for x in hits}
        return next(iter(unique.values())) if len(unique) == 1 else None

    # Fixed-width substitution is the least boundary-ambiguous rescue.
    for hits in (substitution_hits, deletion_hits, insertion_hits):
        hit = unique_hit(hits)
        if hit is not None:
            return hit
    return (None, None, None, 'none')


def iter_CARLIN_raw_read_chunks(data_path, sample, *, cell_barcode_length, umi_length=8, chunk_size=100000):
    """Yield paired DARLIN FASTQ records as bounded DataFrame chunks.

    Each record carries its global ``read_ordinal`` (0-based FASTQ order) so
    downstream stages can index sidecar rows without retaining every read ID
    in memory. Duplicate read-ID detection hashes identifiers to 8-byte
    digests: at Illumera-scale read counts a 64-bit digest collision
    probability stays negligible and a collision only fails closed as a false
    duplicate report, never silently drops a read.
    """
    if cell_barcode_length not in {8, 10}:
        raise ValueError(f'DARLIN cell barcode length must be 8 or 10, found {cell_barcode_length!r}')
    if umi_length != 8:
        raise ValueError(f'DARLIN UMI length must be 8, found {umi_length!r}')
    chunk_size = int(chunk_size)
    if chunk_size < 1:
        raise ValueError(f'DARLIN chunk_size must be >=1, found {chunk_size!r}')
    tag_length = cell_barcode_length + umi_length
    seq_path = f'{data_path}/{sample}_R1.fastq.gz'
    tag_path = f'{data_path}/{sample}_R2.fastq.gz'
    for path in (seq_path, tag_path):
        with open(path, 'rb') as raw_handle:
            if raw_handle.read(2) != b'\x1f\x8b':
                raise ValueError(f'DARLIN FASTQ is not a valid gzip stream: {path}')
    seen_read_id_hashes = set()
    record_number = 0
    batch_lines = 4 * chunk_size

    def _normalize_headers(headers):
        read_ids = []
        for header in headers:
            if not header.startswith('@'):
                raise ValueError(f'Invalid FASTQ header (expected @): {header[:80]!r}')
            read_id = header[1:].split(None, 1)[0].strip()
            if read_id.endswith(('/1', '/2')):
                read_id = read_id[:-2]
            if not read_id:
                raise ValueError('DARLIN read identifier is empty after normalization')
            read_ids.append(read_id)
        return read_ids

    try:
        with gzip.open(seq_path, 'rt', encoding='ascii', newline='') as seq_handle, gzip.open(tag_path, 'rt', encoding='ascii', newline='') as tag_handle:
            seq_iter = iter(seq_handle)
            tag_iter = iter(tag_handle)
            while True:
                seq_lines = list(islice(seq_iter, batch_lines))
                tag_lines = list(islice(tag_iter, batch_lines))
                if not seq_lines and not tag_lines:
                    break
                base = record_number
                if len(seq_lines) != len(tag_lines):
                    record_number += max(len(seq_lines), len(tag_lines)) // 4 + 1
                    raise ValueError(f'DARLIN paired FASTQ record counts differ at record {record_number}: {seq_path} vs {tag_path}')
                if len(seq_lines) % 4:
                    record_number += len(seq_lines) // 4 + 1
                    raise ValueError(f'Incomplete DARLIN FASTQ record {record_number}: {seq_path}')
                count = len(seq_lines) // 4
                record_number += count
                headers = [line.rstrip('\r\n') for line in seq_lines[0::4]]
                sequences = [line.rstrip('\r\n') for line in seq_lines[1::4]]
                pluses = seq_lines[2::4]
                qualities = [line.rstrip('\r\n') for line in seq_lines[3::4]]
                tag_headers = [line.rstrip('\r\n') for line in tag_lines[0::4]]
                tags = [line.rstrip('\r\n') for line in tag_lines[1::4]]
                tag_qualities = [line.rstrip('\r\n') for line in tag_lines[3::4]]
                bad_plus = next((index for index, line in enumerate(pluses) if not line.rstrip('\r\n').startswith('+')), None)
                if bad_plus is not None:
                    raise ValueError(f'Invalid FASTQ + line at record {base + bad_plus + 1}: {seq_path}')
                bad_length = next((index for index, (sequence, quality) in enumerate(zip(sequences, qualities, strict=True)) if len(sequence) != len(quality)), None)
                if bad_length is not None:
                    raise ValueError(
                        f'DARLIN FASTQ sequence/quality length mismatch at record {base + bad_length + 1}: '
                        f'{seq_path} ({len(sequences[bad_length])} != {len(qualities[bad_length])})'
                    )
                seq_ids = _normalize_headers(headers)
                tag_ids = _normalize_headers(tag_headers)
                if seq_ids != tag_ids:
                    index = next(index for index, (left, right) in enumerate(zip(seq_ids, tag_ids, strict=True)) if left != right)
                    raise ValueError(f'DARLIN paired FASTQ read-ID mismatch at record {base + index + 1}: {seq_ids[index]!r} != {tag_ids[index]!r}')
                for read_id, read_hash in zip(seq_ids, (hashlib.blake2b(read_id.encode('ascii'), digest_size=8).digest() for read_id in seq_ids), strict=True):
                    if read_hash in seen_read_id_hashes:
                        raise ValueError(f'Duplicate DARLIN read identifier in sample {sample!r}: {read_id}')
                    seen_read_id_hashes.add(read_hash)
                bad_tag = next((index for index, tag in enumerate(tags) if len(tag) < tag_length), None)
                if bad_tag is not None:
                    raise ValueError(f'DARLIN tag read is shorter than barcode+UMI at record {base + bad_tag + 1}: required {tag_length}, found {len(tags[bad_tag])}')
                cell_bcs = [tag[:cell_barcode_length] for tag in tags]
                umis = [tag[cell_barcode_length:tag_length] for tag in tags]
                # R2 tags share one fixed barcode+UMI width, so the whole batch's
                # Phred summaries collapse into one vectorized array reduction.
                tag_block = ''.join(quality[:tag_length] for quality in tag_qualities)
                tag_phred = np.frombuffer(tag_block.encode('ascii'), dtype=np.uint8).reshape(count, tag_length).astype(np.int32) - 33
                clone_means = []
                clone_mins = []
                for quality in qualities:
                    raw = quality.encode('ascii')
                    if not raw:
                        raise ValueError('Cannot summarize an empty FASTQ quality segment')
                    clone_means.append((sum(raw) - 33 * len(raw)) / len(raw))
                    clone_mins.append(min(raw) - 33)
                yield pd.DataFrame({
                    'read_id': seq_ids,
                    'read_ordinal': range(base, base + count),
                    'cell_bc': cell_bcs,
                    'cell_bc_quality_mean': tag_phred[:, :cell_barcode_length].mean(axis=1).tolist(),
                    'cell_bc_quality_min': tag_phred[:, :cell_barcode_length].min(axis=1).tolist(),
                    'library': [sample] * count,
                    'cell_id': [f'{sample}_{cell_bc}' for cell_bc in cell_bcs],
                    'umi': umis,
                    'umi_quality_mean': tag_phred[:, cell_barcode_length:tag_length].mean(axis=1).tolist(),
                    'umi_quality_min': tag_phred[:, cell_barcode_length:tag_length].min(axis=1).tolist(),
                    'umi_id': [f'{cell_bc}_{umi}' for cell_bc, umi in zip(cell_bcs, umis, strict=True)],
                    'clone_id': sequences,
                    'clone_id_quality_mean': clone_means,
                    'clone_id_quality_min': clone_mins,
                })
    except (gzip.BadGzipFile, EOFError, UnicodeDecodeError) as error:
        raise ValueError(
            f'DARLIN FASTQ gzip stream is corrupt or truncated while reading sample {sample!r} '
            f'around record {record_number}: {type(error).__name__}: {error}'
        ) from error


def _build_unique_hamming1_barcode_index(ref_cell_barcodes):
    """Build fail-closed exact and unique Hamming-1 barcode indexes."""
    canonical = tuple(str(x).strip().upper() for x in ref_cell_barcodes)
    exact = set(canonical)
    unique = {}
    ambiguous = set()
    for barcode in canonical:
        for i, original in enumerate(barcode):
            for base in 'ACGT':
                if base == original:
                    continue
                neighbor = barcode[:i] + base + barcode[i + 1:]
                if neighbor in exact:
                    continue
                previous = unique.get(neighbor)
                if previous is None and neighbor not in ambiguous:
                    unique[neighbor] = barcode
                elif previous != barcode:
                    unique.pop(neighbor, None)
                    ambiguous.add(neighbor)
    return exact, unique, ambiguous


def correct_darlin_cell_barcodes(df_input, ref_cell_barcodes=None, *, barcode_index=None):
    """Apply RNA-compatible unique Hamming-1 correction to DARLIN cell barcodes.

    ``barcode_index`` may be precomputed once per sample to avoid rebuilding the
    Hamming-1 neighborhood for every streaming chunk. Classification is fully
    vectorized: per-read Python loops would dominate chunk-level wall time.
    """
    df = df_input.copy()
    if barcode_index is None:
        if ref_cell_barcodes is None:
            raise ValueError('ref_cell_barcodes or barcode_index is required')
        barcode_index = _build_unique_hamming1_barcode_index(ref_cell_barcodes)
    exact, unique, ambiguous = barcode_index
    observed = df['cell_bc'].astype(str).str.upper()
    is_exact = observed.isin(exact)
    is_ambiguous = ~is_exact & observed.isin(ambiguous)
    corrected_unique = observed.map(unique)
    is_unique = ~is_exact & ~is_ambiguous & corrected_unique.notna()
    status = np.where(is_exact, 'exact', np.where(is_ambiguous, 'ambiguous', np.where(is_unique, 'corrected_hamming1', 'unmatched')))
    corrected = np.where(is_exact, observed.to_numpy(), np.where(is_unique, corrected_unique.astype(str).to_numpy(), ''))
    df['cell_bc_original'] = observed.to_numpy()
    df['cell_bc_correction_status'] = status
    df['cell_bc'] = corrected
    valid = df['cell_bc'].ne('')
    df.loc[valid, 'cell_id'] = df.loc[valid, 'library'].astype(str) + '_' + df.loc[valid, 'cell_bc'].astype(str)
    df.loc[valid, 'umi_id'] = df.loc[valid, 'cell_bc'].astype(str) + '_' + df.loc[valid, 'umi'].astype(str)
    return df


def CARLIN_preprocessing(df_input, template='cCARLIN', ref_cell_barcodes=None, seq_5prime_upper_N=None, seq_3prime_upper_N=None, seq_3prime_max_mismatches=0, seq_5prime_max_edits=1, return_read_level=False, use_short_3prime=None, progress=None):
    """Filter DARLIN reads with indel-aware conserved-anchor matching.

    The editable array is never fuzzy-matched. Only the short conserved 5' and
    3' anchors may absorb sequencing/PCR errors; the extracted allele sequence
    remains the observed molecule sequence and is resolved later at UMI level.
    ``use_short_3prime`` pins the primer-variant decision for callers that
    stream many chunks; ``None`` keeps the historical per-call inference so
    mixed-length streams cannot silently switch primers between chunks.
    ``progress`` defaults to ``print`` for standalone use; streaming callers
    pass a collector so per-chunk diagnostics aggregate instead of printing
    thousands of near-identical blocks.
    """
    if progress is None:
        progress = print
    df_output = df_input.copy()
    tot_fastq_N = len(df_output)
    progress(f'Total fastq: {tot_fastq_N}')
    stable_columns = ['cell_bc', 'library', 'cell_id', 'umi', 'umi_id', 'clone_id', '5prime_match_type', '5prime_edit_distance', '3prime_detected', '3prime_match_type', '3prime_mismatch_count', 'read']
    if tot_fastq_N == 0:
        empty = pd.DataFrame(columns=stable_columns)
        read_level = pd.DataFrame(columns=['read_id', 'cell_bc', 'cell_id', 'umi', 'clone_id', '3prime_detected', '3prime_match_type', '3prime_mismatch_count'])
        return (empty, read_level) if return_read_level else empty

    seq_length = int(df_output['clone_id'].iloc[:100].str.len().mean())
    if use_short_3prime is None:
        use_short_3prime = seq_length < 300
    if use_short_3prime:
        progress('Use short primer')
    seq_5prime, seq_5prime_full, seq_3prime, appendix, CARLIN_seq = _template_parts(template, use_short_3prime=use_short_3prime)
    seq_full = seq_5prime_full + CARLIN_seq + seq_3prime[len(appendix):]
    if seq_5prime_upper_N is not None:
        seq_5prime = seq_5prime[-seq_5prime_upper_N:]
    if seq_3prime_upper_N is None:
        seq_3prime_upper_N = len(seq_3prime)
    seq_3prime_0 = seq_3prime[:seq_3prime_upper_N]

    # Fast exact path first; only misses pay the bounded indel-aware cost.
    df_output['_5prime_start'] = df_output['clone_id'].str.find(seq_5prime)
    df_output['_5prime_end'] = df_output['_5prime_start'] + len(seq_5prime)
    df_output['5prime_edit_distance'] = np.where(df_output['_5prime_start'].ge(0), 0, -1)
    df_output['5prime_match_type'] = np.where(df_output['_5prime_start'].ge(0), 'exact', 'none')
    missing_5 = df_output['_5prime_start'].lt(0)
    if int(seq_5prime_max_edits) > 0 and missing_5.any():
        matches = df_output.loc[missing_5, 'clone_id'].apply(
            lambda read: find_indel_aware_subsequence(read, seq_5prime, max_edits=seq_5prime_max_edits, downstream_reference=CARLIN_seq)
        )
        df_output.loc[missing_5, '_5prime_start'] = matches.map(lambda x: -1 if x[0] is None else int(x[0]))
        df_output.loc[missing_5, '_5prime_end'] = matches.map(lambda x: -1 if x[1] is None else int(x[1]))
        df_output.loc[missing_5, '5prime_edit_distance'] = matches.map(lambda x: -1 if x[2] is None else int(x[2]))
        df_output.loc[missing_5, '5prime_match_type'] = matches.map(lambda x: x[3])
    df_output['Valid_5prime'] = df_output['_5prime_start'].ge(0)

    # Avoid DataFrame.apply(axis=1) on millions of reads.
    df_output['_after_5prime'] = [
        str(seq)[int(end):] if int(end) >= 0 else ''
        for seq, end in zip(df_output['clone_id'].to_numpy(), df_output['_5prime_end'].to_numpy(), strict=True)
    ]
    df_output['_3prime_exact_index'] = df_output['_after_5prime'].str.find(seq_3prime_0)
    df_output['Valid_3prime_0'] = df_output['_3prime_exact_index'].ge(0)
    df_output['Valid_3prime_fuzzy'] = False
    df_output['3prime_fuzzy_index'] = -1
    df_output['3prime_mismatch_count'] = -1
    if seq_3prime_max_mismatches > 0:
        fuzzy_candidates = df_output['Valid_5prime'] & ~df_output['Valid_3prime_0']
        fuzzy_matches = df_output.loc[fuzzy_candidates, '_after_5prime'].apply(
            lambda read: find_approximate_subsequence(read, seq_3prime_0, max_mismatches=seq_3prime_max_mismatches)
        )
        if len(fuzzy_matches) > 0:
            fuzzy_index = fuzzy_matches.map(lambda x: -1 if x[0] is None else x[0])
            fuzzy_mismatch = fuzzy_matches.map(lambda x: -1 if x[1] is None else x[1])
            df_output.loc[fuzzy_candidates, '3prime_fuzzy_index'] = fuzzy_index
            df_output.loc[fuzzy_candidates, '3prime_mismatch_count'] = fuzzy_mismatch
            df_output.loc[fuzzy_candidates, 'Valid_3prime_fuzzy'] = fuzzy_mismatch >= 0

    if seq_length < len(seq_full):
        progress(f'Fastq length insufficient ({seq_length} bp)')
        df_output['Valid_3prime_1'] = df_output['Valid_5prime'] & df_output['clone_id'].str[-seq_3prime_upper_N:].map(seq_full.__contains__)
        df_output['Valid_3prime'] = df_output['Valid_3prime_0'] | df_output['Valid_3prime_fuzzy'] | df_output['Valid_3prime_1']
    else:
        df_output['Valid_3prime_1'] = False
        df_output['Valid_3prime'] = df_output['Valid_3prime_0'] | df_output['Valid_3prime_fuzzy']
    df_output['Valid'] = df_output['Valid_5prime'] & df_output['Valid_3prime']
    progress(f"Fastq frac. with valid 5 prime: {np.mean(df_output['Valid_5prime']):.3f}")
    progress(f"Fastq frac. with valid 3 prime: {np.mean(df_output['Valid_3prime']):.3f}")
    df_output = df_output.query('Valid==True').copy()
    if ref_cell_barcodes is not None:
        df_output = df_output[df_output['cell_bc'].isin(ref_cell_barcodes)]
        valid_BC_N = len(df_output)
        progress(f'Fastq with valid barcodes: {valid_BC_N} ({valid_BC_N / tot_fastq_N:.3f})')

    exact = df_output['Valid_3prime_0']
    fuzzy = df_output['Valid_3prime_fuzzy'] & ~exact
    read_end = df_output['Valid_3prime_1'] & ~exact & ~fuzzy
    selected = exact | fuzzy | read_end
    if not bool(selected.any()):
        df_output = df_output.iloc[0:0].copy()
    else:
        exact, fuzzy, read_end, selected = (series.loc[selected] for series in (exact, fuzzy, read_end, selected))
        df_output = df_output.loc[selected].copy()
        row_count = len(df_output)
        after_values = df_output['_after_5prime'].astype(str).tolist()
        exact_flags = exact.to_numpy()
        fuzzy_flags = fuzzy.to_numpy()
        read_end_flags = read_end.to_numpy()
        exact_indices = df_output['_3prime_exact_index'].to_numpy()
        fuzzy_indices = df_output['3prime_fuzzy_index'].to_numpy()
        fuzzy_mismatches = pd.to_numeric(df_output['3prime_mismatch_count'], errors='coerce').fillna(-1).to_numpy()
        alleles = [None] * row_count
        match_types = [None] * row_count
        mismatch_counts = [-1] * row_count
        detected_flags = [False] * row_count
        for position in range(row_count):
            after = after_values[position]
            if exact_flags[position]:
                alleles[position] = after[:int(exact_indices[position])] + appendix
                match_types[position] = 'exact'
                mismatch_counts[position] = 0
                detected_flags[position] = True
            elif fuzzy_flags[position]:
                alleles[position] = after[:int(fuzzy_indices[position])] + appendix
                match_types[position] = 'fuzzy'
                mismatch_counts[position] = int(fuzzy_mismatches[position])
                detected_flags[position] = True
            elif read_end_flags[position]:
                alleles[position] = after
                match_types[position] = 'read_end'
                mismatch_counts[position] = -1
                detected_flags[position] = False
        df_output['clone_id'] = alleles
        df_output['3prime_detected'] = detected_flags
        df_output['3prime_match_type'] = match_types
        df_output['3prime_mismatch_count'] = mismatch_counts
    progress(f'Fastq frac. with valid 3 and 5 prime: {len(df_output) / tot_fastq_N:.3f}')
    df_output = check_editing(df_output, template, progress=progress)
    read_level_columns = ['read_id', 'cell_bc', 'cell_id', 'umi', 'clone_id', '3prime_detected', '3prime_match_type', '3prime_mismatch_count']
    if 'read_ordinal' in df_output.columns:
        read_level_columns = [*read_level_columns, 'read_ordinal']
    read_level = df_output.reindex(columns=read_level_columns).copy()
    if df_output.empty:
        return (pd.DataFrame(columns=stable_columns), read_level) if return_read_level else pd.DataFrame(columns=stable_columns)
    df_output['unique_id'] = df_output['cell_id'] + '_' + df_output['umi_id'] + '_' + df_output['clone_id']
    read_counts = df_output.groupby('unique_id').agg(read=('unique_id', 'count')).reset_index()
    keep = ['cell_bc', 'library', 'cell_id', 'umi', 'umi_id', 'clone_id', 'unique_id', 'Valid', '5prime_match_type', '5prime_edit_distance', '3prime_detected', '3prime_match_type', '3prime_mismatch_count']
    aggregated = df_output.filter(keep).merge(read_counts, on='unique_id').drop(['Valid', 'unique_id'], axis=1).drop_duplicates()
    return (aggregated, read_level) if return_read_level else aggregated


def denoise_clonal_data(df_raw, read_cutoff=3, group_keys=('library', 'cell_id', 'cell_bc', 'clone_id', 'umi'), return_lineage=False):
    """Denoise clone IDs within biological cells and optionally return lineage."""
    df_input = df_raw.copy()
    df_input['_input_clone_id'] = df_input['clone_id']
    pass_read = df_input['read'] >= read_cutoff
    df_input['_passes_read_cutoff'] = pass_read
    df_input.loc[~pass_read, 'clone_id'] = np.nan
    passing = df_input.loc[pass_read]
    local_keys = [key for key in ('library', 'cell_id', 'cell_bc') if key in df_input.columns]
    if not passing.empty:
        if local_keys:
            for _, index in passing.groupby(local_keys, sort=False, observed=True).groups.items():
                index = list(index)
                df_input.loc[index, 'clone_id'] = denoise_sequence(df_input.loc[index, 'clone_id'], df_input.loc[index, 'read'])
        else:
            df_input.loc[pass_read, 'clone_id'] = denoise_sequence(passing['clone_id'], passing['read'])
    df_input['clone_id'] = df_input['clone_id'].replace('nan', np.nan)
    lineage_keys = [key for key in ('library', 'cell_id', 'cell_bc', 'umi', '3prime_detected', '3prime_match_type', '3prime_mismatch_count') if key in df_input.columns]
    lineage = df_input[lineage_keys + ['_input_clone_id', 'clone_id', '_passes_read_cutoff']].rename(columns={'_input_clone_id': 'input_clone_id', 'clone_id': 'corrected_clone_id', '_passes_read_cutoff': 'passes_read_cutoff'}).drop_duplicates().reset_index(drop=True)
    group_keys = [key for key in group_keys if key in df_input.columns]
    grouped = group_cells(df_input.drop(columns=['_input_clone_id', '_passes_read_cutoff']).dropna(), group_keys=group_keys)
    return (grouped, lineage) if return_lineage else grouped


def denoise_sequence(input_seqs, read_count):
    """Collapse near-identical clone sequences with local directional clustering."""
    seqs = [str(value).encode('ascii') for value in input_seqs]
    counts = list(read_count)
    if not seqs:
        return np.array([], dtype=str)
    if len(counts) != len(seqs):
        raise ValueError('read_count does not have the same size as input_seqs')
    seq_count = {}
    for seq, count in zip(seqs, counts, strict=True):
        seq_count[seq] = seq_count.get(seq, 0) + int(count)
    from umi_tools import UMIClusterer
    mapping = {}
    for cluster in UMIClusterer(cluster_method='directional')(seq_count, threshold=1):
        parent = cluster[0]
        for seq in cluster:
            mapping[seq] = parent
    return np.asarray([mapping[seq].decode('ascii') for seq in seqs], dtype=str)


def add_read_fractions(df_input):
    """Add per-cell and per-clone read fractions used for final filtering."""
    df_output = df_input.drop(columns=['total_read_per_cell', 'total_read_per_clone'], errors='ignore').copy()
    df_output['total_read_per_cell'] = df_output.groupby('cell_id')['read'].transform('max')
    df_output['total_read_per_clone'] = df_output.groupby('clone_id')['read'].transform('max')
    df_output['cell_read_fraction'] = df_output['read'] / df_output['total_read_per_cell']
    df_output['clone_read_fraction'] = df_output['read'] / df_output['total_read_per_clone']
    return df_output


def estimate_read_cutoff(df):
    """Pick the first read cutoff where the cell-count slope stabilizes."""
    df = df.sort_values('read_cutoff').reset_index(drop=True)
    df['difference'] = df['cell_id_count'].diff()
    df['slope_change'] = df['difference'].diff()
    for i in range(1, len(df) - 1):
        current = abs(df.loc[i, 'slope_change'])
        next_val = abs(df.loc[i + 1, 'slope_change'])
        if current < 5 and next_val < 5:
            return df.loc[i, 'read_cutoff']
    return 7


def unique_cell_cutoff_stats(df):
    """Count unique cells and clones across exponentially spaced read cutoffs."""
    base = 1.5
    max_read = df['read'].max() if not df.empty else 1
    upper_log = max(1, int(np.ceil(np.log(max_read) / np.log(base)))) if max_read > 1 else 1
    read_cutoffs = sorted({int(base ** x) for x in range(upper_log)})
    stats = {'read_cutoff': read_cutoffs}
    for key in ('cell_id', 'clone_id'):
        stats[f'{key}_count'] = [df.loc[df['read'] >= read_cutoff, key].nunique() for read_cutoff in read_cutoffs]
    return pd.DataFrame(stats)


def group_cells(df_HQ, group_keys=('library', 'cell_id', 'clone_id')):
    """Aggregate read counts, and UMI counts when UMI is not part of the grouping key."""
    agg_spec = {'read': ('read', 'sum')}
    if 'umi' not in group_keys and 'umi' in df_HQ.columns:
        agg_spec['umi_count'] = ('umi', 'count')
    df_out = df_HQ.groupby(group_keys).agg(**agg_spec)
    return df_out.reset_index()


darlin_yaml = None


ASSIGNMENT_METHOD = 'biological_cell_exact_barcode_union_v3'


REFERENCE_EVIDENCE_COLUMN = 'sample_count'


OBSERVED_EVIDENCE_COLUMN = 'observed_biological_sample_count'


QC_PASSING_EVIDENCE_COLUMN = 'qc_passing_biological_sample_count'


REFERENCE_EVIDENCE_METHOD = 'exact_barcode_qc_passing_distinct_biological_samples_per_locus_v3'


BIOLOGICAL_SOURCE_COLUMN = 'biological_source_sample'


CANONICAL_CELL_BARCODE_COLUMN = 'canonical_cell_barcode'


NORMALIZED_BARCODE_COLUMN = 'normalized_exact_barcode'


SOURCE_DUPLICATE_COUNT_COLUMN = 'source_duplicate_count'


SOURCE_DUPLICATE_WARNING_COLUMN = 'source_duplicate_warning'


SOURCE_DUPLICATE_CONFLICT_COLUMN = 'source_duplicate_conflict'


DEFAULT_TEMPLATE_LOCUS = {'cCARLIN': 'CA', 'Rosa': 'RA', 'Rosa_v2': 'RA', 'Tigre': 'TA', 'Tigre_2022': 'TA', 'Tigre_2022_v2': 'TA'}


FINAL_COLUMNS = ['Sample_ID', 'DARLIN_Clone_ID', 'DARLIN_QC_Pass', 'DARLIN_QC_Status', 'DARLIN_Loci', 'DARLIN_Barcode_Summary', 'DARLIN_Total_Read_Count', 'DARLIN_Total_UMI_Count', 'DARLIN_Min_Barcode_Read_Count', 'DARLIN_Min_Barcode_UMI_Count', 'DARLIN_Min_Cell_Read_Fraction', 'DARLIN_3Prime_Evidence', 'DARLIN_Max_3Prime_Mismatch_Count', 'DARLIN_Max_Ambiguous_Base_Count', 'DARLIN_Barcode_Complexity', 'DARLIN_Assignment_Method', 'DARLIN_Assignment_Read_Support_Fraction', 'DARLIN_Assignment_UMI_Support_Fraction', 'DARLIN_Reference_QC_Status', 'DARLIN_Reference_Evidence_Method', 'DARLIN_Reference_Evidence_Complete', 'DARLIN_Reference_Evidence_Barcodes', 'DARLIN_Max_Reference_Sample_Count_Observed', 'DARLIN_Max_QC_Passing_Biological_Sample_Count', 'DARLIN_Double_Cell_Suspect', 'DARLIN_Double_Cell_Reason', 'DARLIN_Aggregation_Warnings', 'DARLIN_Source_Duplicate_Count', 'DARLIN_Barcodes_Total', 'DARLIN_Barcodes_After_Read_QC', 'DARLIN_Barcodes_After_UMI_QC', 'DARLIN_Barcodes_After_Edit_QC', 'DARLIN_Barcodes_After_3Prime_QC', 'DARLIN_Barcodes_After_3Prime_Mismatch_QC', 'DARLIN_Barcodes_After_Reference_QC', 'DARLIN_Barcodes_After_Sequence_QC', 'DARLIN_Barcodes_Used', 'QC_Extraction_Read_Cutoff', 'QC_Extraction_3Prime_Max_Mismatches', 'QC_Extraction_Min_Cell_Read_Fraction', 'QC_Min_Read_Count', 'QC_Min_UMI_Count', 'QC_Require_Edited_Barcode', 'QC_Require_3Prime_Detected', 'QC_Max_3Prime_Mismatches', 'QC_Allow_Read_End_Rescue', 'QC_Min_Read_End_Read_Count', 'QC_Min_Read_End_UMI_Count', 'QC_Min_Read_End_Cell_Read_Fraction', 'QC_Min_Read_End_Barcode_Length', 'QC_Max_Ambiguous_Bases', 'QC_Max_Barcode_Complexity', 'QC_Max_Reference_Sample_Count', 'QC_Require_Reference_Evidence', 'Barcode_Map_SHA256', 'Barcode_Map_Rows']


BARCODE_COLUMNS = ['library', 'cell_id', 'cell_bc', 'clone_id', '5prime_match_type', '5prime_edit_distance', '3prime_detected', '3prime_match_type', '3prime_mismatch_count', 'read', 'umi_count', 'bc_length', 'total_read_per_cell', 'total_read_per_clone', 'cell_read_fraction', 'clone_read_fraction', 'edited', REFERENCE_EVIDENCE_COLUMN]


AGGREGATION_COLUMNS = ['evidence_class', 'source_sample', BIOLOGICAL_SOURCE_COLUMN, CANONICAL_CELL_BARCODE_COLUMN, NORMALIZED_BARCODE_COLUMN, 'source_samples', 'physical_source_count', SOURCE_DUPLICATE_COUNT_COLUMN, SOURCE_DUPLICATE_WARNING_COLUMN, SOURCE_DUPLICATE_CONFLICT_COLUMN, OBSERVED_EVIDENCE_COLUMN, QC_PASSING_EVIDENCE_COLUMN, 'template', 'locus']


def _load_runtime_deps():
    """Import heavy table dependencies only when work is executed."""
    global pd, darlin_yaml
    if darlin_yaml is None:
        import yaml as _yaml
        darlin_yaml = _yaml


def _read_yaml(path):
    _load_runtime_deps()
    if path is None:
        return {}
    if not Path(path).exists():
        return {}
    with open(path, encoding='utf-8') as handle:
        return select_config(darlin_yaml.safe_load(handle) or {})


def _clone_config(config):
    section = config.get('DARLIN') or {}
    configured_qc = section.get('qc') or {}
    unknown_qc = sorted(set(configured_qc).difference(DEFAULT_QC))
    if unknown_qc:
        raise ValueError('DARLIN.qc contains unknown key(s): ' + ', '.join(unknown_qc))
    qc = dict(DEFAULT_QC)
    qc.update(configured_qc)
    template_locus = dict(DEFAULT_TEMPLATE_LOCUS)
    template_locus.update(section.get('template_locus') or {})
    extraction = dict(DEFAULT_EXTRACTION)
    extraction.update(section.get('extraction') or {})
    return {'template_locus': template_locus, 'extraction': extraction, 'qc': qc}


def _ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def _read_csv_rows(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        delimiter = '\t' if Path(path).suffix.lower() in {'.tsv', '.txt'} else ','
        return [{k: (v or '').strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def _resolve_project_path(value, project_dir):
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else project_dir / path


def _project_dir(raw_config):
    config_path = Path(str(raw_config.get('__configfile__', ''))).expanduser()
    if config_path.name == 'config.yaml' and config_path.parent.name == '00_config':
        return config_path.resolve().parent.parent
    return Path.cwd().resolve()


def _cell_name_maps(raw_config):
    project_dir = _project_dir(raw_config)
    manifest = _resolve_project_path(raw_config.get('sample_manifest', '00_config/sample_manifest.tsv'), project_dir)
    barcode_map = _resolve_project_path(raw_config.get('barcode_map', '00_config/Barcode_Map.csv'), project_dir)
    for label, path in (('项目样本清单', manifest), ('项目条码表', barcode_map)):
        if not Path(path).is_file():
            raise FileNotFoundError(f'{label}不存在: {path}')
    project_by_sample = {}
    for row in _read_csv_rows(manifest):
        project_sample = row.get('sample_id', '')
        if not project_sample:
            continue
        for key in ('sample_id', 'rna_sample', 'darlin_sample'):
            value = row.get(key, '')
            if value:
                project_by_sample[value] = project_sample
    plate_by_rna_barcode = {row.get('RNA_Barcode', ''): row.get('PlateID', '') for row in _read_csv_rows(barcode_map) if row.get('RNA_Barcode', '') and row.get('PlateID', '')}
    return (project_by_sample, plate_by_rna_barcode)


def _biological_source_sample(source_sample):
    """Collapse physical RNA-derived inputs onto their biological library."""
    return str(source_sample).strip().removesuffix('__RNA')


def _canonical_cell_id(cell_id, *, source_sample, cell_barcode, project_by_sample, plate_by_rna_barcode):
    """Return one biological cell identity across raw and RNA-derived inputs."""
    physical_source = str(source_sample).strip()
    biological_source = _biological_source_sample(physical_source)
    cell_bc = str(cell_barcode).strip()
    raw_cell_id = str(cell_id).strip()
    if not cell_bc and physical_source and raw_cell_id.startswith(f'{physical_source}_'):
        cell_bc = raw_cell_id[len(physical_source) + 1:]
    project_sample = project_by_sample.get(biological_source, '')
    if not project_sample:
        for candidate in sorted(set(project_by_sample.values()), key=len, reverse=True):
            if candidate and biological_source.startswith(f'{candidate}_'):
                project_sample = candidate
                break
    if not project_sample:
        project_sample = biological_source
    plate_id = plate_by_rna_barcode.get(cell_bc, '')
    if project_sample and plate_id:
        return f'{project_sample}_{plate_id}'
    if project_sample and cell_bc:
        return f'{project_sample}_{cell_bc}'
    if physical_source.endswith('__RNA') and raw_cell_id.startswith(f'{physical_source}_'):
        return f'{biological_source}_{raw_cell_id[len(physical_source) + 1:]}'
    return raw_cell_id


def _collapse_duplicate_sample_ids(final):
    if final.empty or not final['Sample_ID'].duplicated().any():
        return final
    duplicates = sorted(final.loc[final['Sample_ID'].duplicated(keep=False), 'Sample_ID'].astype(str).unique())
    raise ValueError('DARLIN biological cell canonicalization produced duplicate Sample_ID values; refusing to discard evidence: ' + ', '.join(duplicates[:10]))


def _as_bool_series(series):
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({'1', 'true', 't', 'yes', 'y'})


def _numeric_series(df, column, default=0):
    if column not in df.columns:
        return pd.Series(default, index=df.index, dtype='float64')
    return pd.to_numeric(df[column], errors='coerce').fillna(default)


def _read_table(path):
    """读取当前版本的带表头 barcode CSV；损坏或零字节输入直接失败。"""
    _load_runtime_deps()
    df = pd.read_csv(path)
    unnamed = [col for col in df.columns if str(col).startswith('Unnamed:')]
    if unnamed:
        df = df.drop(columns=unnamed)
    return df


def _normalized_exact_barcode(value):
    """Normalize presentation only; never fuzzy-match distinct edit sequences."""
    text = '' if pd.isna(value) else str(value).strip()
    return text.upper()


def _canonical_cell_barcode(row, source_sample):
    value = row.get('cell_bc', '')
    cell_barcode = '' if pd.isna(value) else str(value).strip().upper()
    if cell_barcode:
        return cell_barcode
    raw_cell_id = row.get('cell_id', '')
    raw_cell_id = '' if pd.isna(raw_cell_id) else str(raw_cell_id).strip()
    physical_source = str(source_sample).strip()
    for prefix in (physical_source, _biological_source_sample(physical_source)):
        if prefix and raw_cell_id.startswith(f'{prefix}_'):
            return raw_cell_id[len(prefix) + 1:].strip().upper()
    return raw_cell_id.upper()


def _stable_row_signature(row, columns):
    return tuple(('' if pd.isna(row.get(column, '')) else str(row.get(column, '')).strip() for column in columns))


def _catalog_match_rna_allele(candidate, catalog, *, max_edits=1, min_partial_length=40):
    """Return a unique targeted allele compatible with one RNA-derived call."""
    candidate = _normalized_exact_barcode(candidate)
    if not candidate or not catalog:
        return (None, 'none')
    exact = [target for target in catalog if target == candidate]
    if len(exact) == 1:
        return (exact[0], 'RNA_KNOWN_EXACT')
    compatible = []
    for target in catalog:
        if len(candidate) >= min_partial_length and len(target) >= min_partial_length:
            if target.startswith(candidate) or candidate.startswith(target):
                compatible.append((target, 'RNA_KNOWN_PARTIAL'))
                continue
        if abs(len(candidate) - len(target)) <= max_edits:
            distance = _bounded_edit_distance(candidate, target, max_edits)
            if distance <= max_edits:
                compatible.append((target, 'RNA_KNOWN_EDIT_RESCUE'))
    targets = sorted({target for target, _ in compatible})
    if len(targets) != 1:
        return (None, 'ambiguous' if targets else 'none')
    target = targets[0]
    kinds = {kind for value, kind in compatible if value == target}
    priority = ('RNA_KNOWN_PARTIAL', 'RNA_KNOWN_EDIT_RESCUE')
    kind = next((name for name in priority if name in kinds), sorted(kinds)[0])
    return (target, kind)


def _apply_rna_allele_catalog_rescue(merged, extraction):
    """Gate sparse RNA-derived calls against the paired targeted allele catalog.

    Targeted calls define a sample/locus catalog. RNA-derived evidence may be
    accepted at low read support only when it maps uniquely to that catalog.
    RNA-only novel calls require independent multi-read and multi-UMI support.
    """
    if merged.empty or 'source_sample' not in merged.columns:
        return merged
    merged = merged.copy()
    merged['evidence_class'] = np.where(
        merged['source_sample'].astype(str).str.endswith('__RNA'),
        'RNA_UNCLASSIFIED', 'TARGETED'
    )
    targeted = merged[~merged['source_sample'].astype(str).str.endswith('__RNA')]
    rna_mask = merged['source_sample'].astype(str).str.endswith('__RNA')
    if not rna_mask.any():
        return merged
    max_edits = int(extraction.get('rna_catalog_max_edits', 1))
    min_partial = int(extraction.get('rna_catalog_min_partial_length', 40))
    novel_min_reads = int(extraction.get('rna_novel_min_reads', 2))
    novel_min_umis = int(extraction.get('rna_novel_min_umis', 2))
    # Build the per-(biological sample, locus) targeted allele catalog once
    # instead of re-filtering the whole targeted table inside every row loop.
    catalogs = {
        (biological, locus): sorted(set(group[NORMALIZED_BARCODE_COLUMN].astype(str)))
        for (biological, locus), group in targeted.groupby([BIOLOGICAL_SOURCE_COLUMN, 'locus'], sort=False, observed=True)
    }
    read_support_values = pd.to_numeric(merged['read'], errors='coerce').fillna(0.0) if 'read' in merged.columns else pd.Series(0.0, index=merged.index)
    umi_support_values = pd.to_numeric(merged['umi_count'], errors='coerce').fillna(0.0) if 'umi_count' in merged.columns else pd.Series(0.0, index=merged.index)
    fiveprime_edit_values = pd.to_numeric(merged['5prime_edit_distance'], errors='coerce').fillna(-1.0) if '5prime_edit_distance' in merged.columns else pd.Series(-1.0, index=merged.index)
    detected_values = merged['3prime_detected'].fillna(False).astype(bool) if '3prime_detected' in merged.columns else pd.Series(False, index=merged.index)
    for index in merged.index[rna_mask]:
        row = merged.loc[index]
        biological = row[BIOLOGICAL_SOURCE_COLUMN]
        locus = row['locus']
        catalog = catalogs.get((biological, locus), [])
        target, evidence = _catalog_match_rna_allele(
            row[NORMALIZED_BARCODE_COLUMN], catalog,
            max_edits=max_edits, min_partial_length=min_partial,
        )
        if target is not None:
            merged.at[index, NORMALIZED_BARCODE_COLUMN] = target
            merged.at[index, 'clone_id'] = target
            merged.at[index, 'evidence_class'] = evidence
            continue
        read_support = float(read_support_values.loc[index])
        umi_support = float(umi_support_values.loc[index])
        dual_anchor = bool(detected_values.loc[index]) and float(fiveprime_edit_values.loc[index]) >= 0
        if read_support >= novel_min_reads and umi_support >= novel_min_umis and dual_anchor:
            merged.at[index, 'evidence_class'] = 'RNA_NOVEL_HIGH_CONFIDENCE'
        else:
            merged.at[index, 'evidence_class'] = 'RNA_LOW_EVIDENCE'
    return merged


def _collapse_source_candidates(merged, payload_columns):
    """Deduplicate one physical source without inflating its evidence.

    Exact duplicate rows are represented once. Multiple distinct evidence rows
    for the same physical-source/cell/locus/exact-barcode key are retained as a
    fail-closed marker and represented by their deterministic maximum read/UMI
    values; downstream assignment emits ``SOURCE_DUPLICATE_CONFLICT``.
    """
    physical_key = ['source_sample', BIOLOGICAL_SOURCE_COLUMN, CANONICAL_CELL_BARCODE_COLUMN, 'locus', NORMALIZED_BARCODE_COLUMN]
    rows = []
    for _, group in merged.groupby(physical_key, sort=True, dropna=False):
        signatures = group.apply(lambda row: _stable_row_signature(row, payload_columns), axis=1)
        unique_indices = signatures.drop_duplicates(keep='first').index
        unique = group.loc[unique_indices].copy()
        duplicate_count = int(len(group) - len(unique))
        ordered_indices = sorted(unique.index, key=lambda index: _stable_row_signature(unique.loc[index], payload_columns))
        representative = unique.loc[ordered_indices[0]].copy()
        representative['read'] = float(_numeric_series(unique, 'read').max())
        representative['umi_count'] = float(_numeric_series(unique, 'umi_count').max())
        representative[SOURCE_DUPLICATE_COUNT_COLUMN] = duplicate_count
        representative[SOURCE_DUPLICATE_WARNING_COLUMN] = duplicate_count > 0
        representative[SOURCE_DUPLICATE_CONFLICT_COLUMN] = len(unique) > 1
        rows.append(representative)
    return pd.DataFrame(rows) if rows else merged.iloc[0:0].copy()


def _aggregate_biological_candidates(source_candidates):
    """Sum independent raw/RNA evidence for the same exact biological call."""
    biological_key = [BIOLOGICAL_SOURCE_COLUMN, CANONICAL_CELL_BARCODE_COLUMN, 'locus', NORMALIZED_BARCODE_COLUMN]
    rows = []
    for _, group in source_candidates.groupby(biological_key, sort=True, dropna=False):
        payload_columns = sorted((str(column) for column in group.columns))
        representative_index = min(group.index, key=lambda index: _stable_row_signature(group.loc[index], payload_columns))
        representative = group.loc[representative_index].copy()
        source_samples = sorted({str(value).strip() for value in group['source_sample'] if str(value).strip()})
        biological_source = str(representative[BIOLOGICAL_SOURCE_COLUMN]).strip()
        cell_barcode = str(representative[CANONICAL_CELL_BARCODE_COLUMN]).strip()
        exact_barcode = str(representative[NORMALIZED_BARCODE_COLUMN]).strip()
        representative['source_sample'] = biological_source
        representative['source_samples'] = ';'.join(source_samples)
        representative['physical_source_count'] = len(source_samples)
        representative['library'] = biological_source
        representative['cell_bc'] = cell_barcode
        representative['cell_id'] = f'{biological_source}_{cell_barcode}' if biological_source and cell_barcode else str(representative.get('cell_id', '')).strip()
        representative['clone_id'] = exact_barcode
        representative['read'] = float(_numeric_series(group, 'read').sum())
        representative['umi_count'] = float(_numeric_series(group, 'umi_count').sum())
        representative['bc_length'] = len(exact_barcode)
        representative['edited'] = bool(_as_bool_series(group['edited']).any() if 'edited' in group.columns else False)
        representative['3prime_detected'] = bool(_as_bool_series(group['3prime_detected']).any() if '3prime_detected' in group.columns else False)
        representative['3prime_mismatch_count'] = float(_numeric_series(group, '3prime_mismatch_count').min())
        match_types = sorted({str(value).strip() for value in group.get('3prime_match_type', pd.Series('', index=group.index)) if str(value).strip()})
        representative['3prime_match_type'] = '+'.join(match_types) or 'unknown'
        evidence_classes = sorted({str(value).strip() for value in group.get('evidence_class', pd.Series('', index=group.index)) if str(value).strip()})
        representative['evidence_class'] = '+'.join(evidence_classes) or 'TARGETED'
        representative[SOURCE_DUPLICATE_COUNT_COLUMN] = int(_numeric_series(group, SOURCE_DUPLICATE_COUNT_COLUMN).sum())
        representative[SOURCE_DUPLICATE_WARNING_COLUMN] = bool(_as_bool_series(group[SOURCE_DUPLICATE_WARNING_COLUMN]).any())
        representative[SOURCE_DUPLICATE_CONFLICT_COLUMN] = bool(_as_bool_series(group[SOURCE_DUPLICATE_CONFLICT_COLUMN]).any())
        rows.append(representative)
    aggregated = pd.DataFrame(rows) if rows else source_candidates.iloc[0:0].copy()
    if aggregated.empty:
        return aggregated
    cell_group = [BIOLOGICAL_SOURCE_COLUMN, CANONICAL_CELL_BARCODE_COLUMN, 'locus']
    aggregated['total_read_per_cell'] = aggregated.groupby(cell_group, sort=False)['read'].transform('sum')
    aggregated['total_read_per_clone'] = aggregated.groupby(['locus', NORMALIZED_BARCODE_COLUMN], sort=False)['read'].transform('sum')
    aggregated['cell_read_fraction'] = aggregated['read'] / aggregated['total_read_per_cell']
    aggregated['clone_read_fraction'] = aggregated['read'] / aggregated['total_read_per_clone']
    return aggregated


def _attach_recurrence_evidence(aggregated, qc):
    valid = aggregated[NORMALIZED_BARCODE_COLUMN].astype('string').str.strip().ne('')
    recurrence_key = ['locus', NORMALIZED_BARCODE_COLUMN]
    observed = aggregated.loc[valid].groupby(recurrence_key, sort=False)[BIOLOGICAL_SOURCE_COLUMN].nunique(dropna=True).to_dict()
    non_recurrence_pass = _non_recurrence_qc_mask(aggregated, qc)
    if SOURCE_DUPLICATE_CONFLICT_COLUMN in aggregated.columns:
        non_recurrence_pass &= ~_as_bool_series(aggregated[SOURCE_DUPLICATE_CONFLICT_COLUMN])
    conflicting_cell_locus = aggregated.groupby([BIOLOGICAL_SOURCE_COLUMN, CANONICAL_CELL_BARCODE_COLUMN, 'locus'], sort=False)[NORMALIZED_BARCODE_COLUMN].transform('nunique').gt(1)
    non_recurrence_pass &= ~conflicting_cell_locus
    passing = aggregated.loc[valid & non_recurrence_pass].groupby(recurrence_key, sort=False)[BIOLOGICAL_SOURCE_COLUMN].nunique(dropna=True).to_dict()
    keys = list(zip(aggregated['locus'].astype(str), aggregated[NORMALIZED_BARCODE_COLUMN].astype(str), strict=True))
    aggregated[OBSERVED_EVIDENCE_COLUMN] = pd.array([observed.get(key, pd.NA) for key in keys], dtype='Int64')
    aggregated[QC_PASSING_EVIDENCE_COLUMN] = pd.array([passing.get(key, 0) if key in observed else pd.NA for key in keys], dtype='Int64')
    aggregated[REFERENCE_EVIDENCE_COLUMN] = aggregated[QC_PASSING_EVIDENCE_COLUMN]
    return aggregated


def _write_merged_table(frame, output, checksum_output=None):
    destination = Path(output)
    _ensure_parent(destination)
    if destination.suffix == '.gz':
        with destination.open('wb') as raw_handle:
            with gzip.GzipFile(filename='', fileobj=raw_handle, mode='wb', mtime=0) as gzip_handle:
                with io.TextIOWrapper(gzip_handle, encoding='utf-8', newline='') as text_handle:
                    frame.to_csv(text_handle, index=False)
    else:
        frame.to_csv(destination, index=False)
    if checksum_output is not None:
        checksum_path = Path(checksum_output)
        _ensure_parent(checksum_path)
        digest = hashlib.sha256(destination.read_bytes()).hexdigest()
        checksum_path.write_text(f'{digest}  {destination.name}\n', encoding='ascii')


def merge_barcodes(input_files, samples, template, locus, output, *, configfile=None, checksum_output=None):
    """Merge one locus and attach auditable cross-sample recurrence evidence.

    Evidence is first deduplicated within one physical source, then raw DARLIN
    and RNA-derived sources are summed by biological sample, cell barcode,
    locus and exact barcode. ``sample_count`` remains a compatibility alias for
    ``qc_passing_biological_sample_count``; the ungated observed count is kept
    separately for audit.
    """
    _load_runtime_deps()
    if len(input_files) != len(samples):
        raise ValueError('--input-files and --samples must have the same length')
    clone_config = _clone_config(_read_yaml(configfile))
    qc = clone_config['qc']
    extraction = clone_config['extraction']
    frames = []
    for sample, path in zip(samples, input_files, strict=True):
        df = _read_table(path)
        payload_columns = [column for column in BARCODE_COLUMNS if column in df.columns]
        required = {'cell_id', 'clone_id', 'read', 'umi_count'}
        missing = sorted(required.difference(df.columns))
        if missing:
            raise ValueError(f'{path} is missing required column(s): {missing}')
        df['source_sample'] = sample
        df[BIOLOGICAL_SOURCE_COLUMN] = _biological_source_sample(sample)
        df['template'] = template
        df['locus'] = locus
        df[CANONICAL_CELL_BARCODE_COLUMN] = df.apply(lambda row, current_sample=sample: _canonical_cell_barcode(row, current_sample), axis=1)
        df[NORMALIZED_BARCODE_COLUMN] = df['clone_id'].map(_normalized_exact_barcode)
        df.attrs['payload_columns'] = payload_columns
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=BARCODE_COLUMNS + AGGREGATION_COLUMNS)
    if merged.empty:
        for column in (OBSERVED_EVIDENCE_COLUMN, QC_PASSING_EVIDENCE_COLUMN, REFERENCE_EVIDENCE_COLUMN):
            merged[column] = pd.Series(dtype='Int64')
    else:
        merged = _apply_rna_allele_catalog_rescue(merged, extraction)
        payload_columns = [column for column in BARCODE_COLUMNS if column in merged.columns]
        source_candidates = _collapse_source_candidates(merged, payload_columns)
        merged = _aggregate_biological_candidates(source_candidates)
        merged = _attach_recurrence_evidence(merged, qc)
        merged = merged.sort_values([BIOLOGICAL_SOURCE_COLUMN, CANONICAL_CELL_BARCODE_COLUMN, 'locus', NORMALIZED_BARCODE_COLUMN, 'source_samples'], kind='mergesort').reset_index(drop=True)
    _write_merged_table(merged, output, checksum_output)


def _barcode_value(row):
    return str(row.get('clone_id', '')).strip()


def _reference_qc_masks(df, qc):
    required = bool(qc['require_reference_evidence'])
    if REFERENCE_EVIDENCE_COLUMN in df.columns:
        sample_count = pd.to_numeric(df[REFERENCE_EVIDENCE_COLUMN], errors='coerce')
    else:
        sample_count = pd.Series(float('nan'), index=df.index, dtype='float64')
    evidence_available = sample_count.notna() & sample_count.ge(1)
    within_limit = evidence_available & sample_count.le(float(qc['max_reference_sample_count']))
    reference_mask = within_limit if required else within_limit | ~evidence_available
    return (reference_mask, evidence_available)


def _non_recurrence_qc_masks(df, qc):
    read_mask = _numeric_series(df, 'read') >= float(qc['min_read_count'])
    umi_mask = _numeric_series(df, 'umi_count') >= float(qc['min_umi_count'])
    if qc['require_edited_barcode']:
        edit_mask = _as_bool_series(df['edited']) if 'edited' in df.columns else pd.Series(False, index=df.index)
    else:
        edit_mask = pd.Series(True, index=df.index)
    if qc['require_3prime_detected']:
        prime_col = '3prime_detected'
        prime_mask = _as_bool_series(df[prime_col]) if prime_col in df.columns else pd.Series(False, index=df.index)
        if qc.get('allow_read_end_rescue', False):
            match_type = df['3prime_match_type'].astype(str) if '3prime_match_type' in df.columns else pd.Series('', index=df.index)
            barcode_length = df.apply(_barcode_value, axis=1).str.len()
            read_end_mask = match_type.eq('read_end') & (_numeric_series(df, 'read') >= float(qc['min_read_end_read_count'])) & (_numeric_series(df, 'umi_count') >= float(qc['min_read_end_umi_count'])) & (_numeric_series(df, 'cell_read_fraction') >= float(qc['min_read_end_cell_read_fraction'])) & (barcode_length >= int(qc['min_read_end_barcode_length']))
            prime_mask = prime_mask | read_end_mask
    else:
        prime_mask = pd.Series(True, index=df.index)
    mismatch_mask = _numeric_series(df, '3prime_mismatch_count', default=0) <= float(qc['max_3prime_mismatches'])
    max_ambiguous_bases = int(qc['max_ambiguous_bases'])
    ambiguous_count = df.apply(_barcode_value, axis=1).astype('string').str.upper().str.count('N')
    sequence_mask = ambiguous_count <= max_ambiguous_bases
    return (read_mask, umi_mask, edit_mask, prime_mask, mismatch_mask, sequence_mask)


def _non_recurrence_qc_mask(df, qc):
    masks = _non_recurrence_qc_masks(df, qc)
    combined = pd.Series(True, index=df.index)
    for mask in masks:
        combined &= mask
    return combined


def _stage_masks(df, qc):
    read_mask, umi_mask, edit_mask, prime_mask, mismatch_mask, sequence_mask = _non_recurrence_qc_masks(df, qc)
    reference_mask, reference_evidence = _reference_qc_masks(df, qc)
    return (read_mask, umi_mask, edit_mask, prime_mask, mismatch_mask, reference_mask, sequence_mask, reference_evidence)


def _select_barcodes(df, locus_order):
    if df.empty:
        return df
    selected = df.copy()
    selected['DARLIN_Barcode_Value'] = selected.apply(_barcode_value, axis=1)
    selected['DARLIN_Read_Count'] = _numeric_series(selected, 'read')
    selected['DARLIN_UMI_Count'] = _numeric_series(selected, 'umi_count')
    selected['DARLIN_Cell_Read_Fraction'] = _numeric_series(selected, 'cell_read_fraction')
    selected['DARLIN_3Prime_Mismatch_Count'] = _numeric_series(selected, '3prime_mismatch_count')
    if '3prime_match_type' in selected.columns:
        selected['DARLIN_3Prime_Match_Type'] = selected['3prime_match_type'].astype(str)
    else:
        selected['DARLIN_3Prime_Match_Type'] = 'unknown'
    selected['DARLIN_Locus_Order'] = selected['locus'].map(locus_order).fillna(len(locus_order))
    selected = selected.sort_values(['DARLIN_Locus_Order', 'DARLIN_Read_Count', 'DARLIN_UMI_Count', 'DARLIN_Cell_Read_Fraction', 'DARLIN_3Prime_Mismatch_Count', 'DARLIN_Barcode_Value', 'DARLIN_3Prime_Match_Type'], ascending=[True, False, False, False, True, True, True], kind='mergesort')
    return selected.groupby('locus', as_index=False, sort=False).head(1)


def _conflicting_loci(df, locus_order):
    """Identify loci with more than one observed exact barcode for one cell."""
    if df.empty:
        return ()
    conflicts = []
    for locus, locus_group in df.groupby('locus', sort=True):
        values = {str(value).strip() for value in locus_group.apply(_barcode_value, axis=1) if str(value).strip()}
        if len(values) > 1:
            conflicts.append(str(locus))
    return tuple(conflicts)


def _clone_id(signature):
    digest = hashlib.sha1(signature.encode('utf-8')).hexdigest()[:12]
    return f'DARLIN_Clone_{digest}'


def _summarize_selected(selected, ordered_loci):
    if selected.empty:
        return ('', '', '')
    selected = selected.copy()
    selected['DARLIN_Locus_Order'] = selected['locus'].map(ordered_loci).fillna(len(ordered_loci))
    selected = selected.sort_values(['DARLIN_Locus_Order', 'locus', 'DARLIN_Barcode_Value'])
    loci = ','.join(selected['locus'].astype(str).tolist())
    signature_parts = [f'{row.locus}:{row.DARLIN_Barcode_Value}' for row in selected.itertuples(index=False)]
    summary_parts = [f'{row.locus}={row.DARLIN_Barcode_Value}|read={int(row.DARLIN_Read_Count)}|umi={int(row.DARLIN_UMI_Count)}|3prime={row.DARLIN_3Prime_Match_Type}:{int(row.DARLIN_3Prime_Mismatch_Count)}' for row in selected.itertuples(index=False)]
    return (loci, ';'.join(signature_parts), ';'.join(summary_parts))


def _ambiguous_base_count(series):
    if series.empty:
        return 0
    return int(series.astype(str).str.upper().str.count('N').max())


def _summarize_3prime_evidence(selected):
    if selected.empty:
        return ('', 0)
    evidence = [f'{row.locus}:{row.DARLIN_3Prime_Match_Type}:{int(row.DARLIN_3Prime_Mismatch_Count)}' for row in selected.itertuples(index=False)]
    return (';'.join(evidence), int(selected['DARLIN_3Prime_Mismatch_Count'].max()))


def _support_fraction(selected, eligible, column):
    denominator = float(_numeric_series(eligible, column).sum())
    if denominator <= 0:
        return 0.0
    numerator = float(_numeric_series(selected, column).sum())
    return numerator / denominator


def _max_reference_sample_count(candidates, column=REFERENCE_EVIDENCE_COLUMN):
    if candidates.empty or column not in candidates.columns:
        return ''
    values = pd.to_numeric(candidates[column], errors='coerce').dropna()
    if values.empty:
        return ''
    maximum = float(values.max())
    return int(maximum) if maximum.is_integer() else maximum


def _reference_qc_status(candidates_n, evidence_n, passing_n, require_reference_evidence):
    if candidates_n == 0:
        return 'NOT_EVALUATED'
    if evidence_n == 0:
        return 'MISSING_REQUIRED_EVIDENCE' if require_reference_evidence else 'NOT_REQUIRED'
    if passing_n == 0:
        return 'FAIL'
    if passing_n < candidates_n:
        return 'PASS_FILTERED'
    return 'PASS'


def _status(total, read_n, umi_n, edit_n, prime_n, mismatch_n, reference_evidence_n, reference_n, sequence_n, used_n, complexity, qc):
    if total == 0:
        return 'NO_BARCODE'
    if read_n == 0:
        return 'NO_BARCODE_PASSING_READ_QC'
    if umi_n == 0:
        return 'NO_BARCODE_PASSING_UMI_QC'
    if qc['require_edited_barcode'] and edit_n == 0:
        return 'NO_EDITED_BARCODE'
    if qc['require_3prime_detected'] and prime_n == 0:
        return 'NO_3PRIME_BARCODE'
    if mismatch_n == 0:
        return 'TOO_MANY_3PRIME_MISMATCHES'
    if qc['require_reference_evidence'] and reference_evidence_n == 0:
        return 'MISSING_REFERENCE_EVIDENCE'
    if reference_n == 0:
        return 'NO_BARCODE_PASSING_REFERENCE_QC'
    if sequence_n == 0:
        return 'AMBIGUOUS_BARCODE'
    if complexity >= int(qc['max_barcode_complexity']):
        return 'TOO_MANY_BARCODES'
    if used_n == 0:
        return 'NO_BARCODE_USED'
    return 'PASS'


def _verify_checksum_sidecar(path, checksum_path):
    fields = Path(checksum_path).read_text(encoding='ascii').strip().split()
    if not fields:
        raise ValueError(f'empty DARLIN merged-table checksum: {checksum_path}')
    expected = fields[0].lower()
    observed = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if expected != observed:
        raise ValueError(f'DARLIN merged-table checksum mismatch for {path}: expected {expected}, observed {observed}')


def build_clone_table(barcode_files, templates, loci, configfile, output, *, checksum_files=None):
    """Integrate template-level barcodes into the final DARLIN clone table."""
    _load_runtime_deps()
    if len(barcode_files) != len(templates) or len(templates) != len(loci):
        raise ValueError('barcode_files, templates, and loci must have the same length')
    if checksum_files is not None:
        if len(checksum_files) != len(barcode_files):
            raise ValueError('checksum_files and barcode_files must have the same length')
        for path, checksum_path in zip(barcode_files, checksum_files, strict=True):
            _verify_checksum_sidecar(path, checksum_path)
    raw_config = _read_yaml(configfile)
    raw_config['__configfile__'] = str(Path(configfile).expanduser().resolve())
    config = _clone_config(raw_config)
    project_by_sample, plate_by_rna_barcode = _cell_name_maps(raw_config)
    barcode_map_path = _resolve_project_path(raw_config.get('barcode_map', '00_config/Barcode_Map.csv'), _project_dir(raw_config))
    barcode_map_sha256, barcode_map_rows = _barcode_map_digest(barcode_map_path)
    qc = config['qc']
    extraction = config['extraction']
    extraction_read_cutoff = extraction.get('read_cutoff', '')
    extraction_3prime_max_mismatches = extraction.get('seq_3prime_max_mismatches', '')
    extraction_min_cell_read_fraction = extraction.get('min_cell_read_fraction', '')
    frames = []
    for path, template, locus in zip(barcode_files, templates, loci, strict=True):
        df = _read_table(path)
        if df.empty:
            continue
        if 'cell_id' not in df.columns:
            raise ValueError(f'{path} is missing required column: cell_id')
        df['template'] = template
        df['locus'] = locus
        frames.append(df)
    if frames:
        all_barcodes = pd.concat(frames, ignore_index=True)
    else:
        all_barcodes = pd.DataFrame(columns=['cell_id', 'locus', 'clone_id'])
    locus_order = {locus: index for index, locus in enumerate(loci)}
    if all_barcodes.empty:
        all_barcodes['_canonical_cell_id'] = pd.Series(dtype='string')
    else:
        # Fail closed before canonicalization: a cell barcode absent from
        # Barcode_Map.csv must never silently fall back to a <sample>_<barcode>
        # identity that diverges from RNA.h5ad obs_names.
        observed_cell_barcodes = {str(barcode) for barcode in all_barcodes['cell_bc'].dropna().astype(str) if str(barcode).strip()}
        unmapped_barcodes = sorted(observed_cell_barcodes.difference(plate_by_rna_barcode))
        if unmapped_barcodes:
            raise ValueError(
                'DARLIN cell barcodes absent from Barcode_Map.csv; refusing to emit '
                'Sample_ID values divergent from RNA cell identities: '
                + ', '.join(unmapped_barcodes[:10])
            )
        all_barcodes['_canonical_cell_id'] = all_barcodes.apply(lambda row: _canonical_cell_id(row.get('cell_id', ''), source_sample=row.get('source_sample', ''), cell_barcode=row.get('cell_bc', ''), project_by_sample=project_by_sample, plate_by_rna_barcode=plate_by_rna_barcode), axis=1)
    rows = []
    for standard_sample_id, group in all_barcodes.groupby('_canonical_cell_id', sort=True):
        group = group.copy()
        read_mask, umi_mask, edit_mask, prime_mask, mismatch_mask, reference_mask, sequence_mask, reference_evidence = _stage_masks(group, qc)
        read_pass = group[read_mask]
        umi_pass = group[read_mask & umi_mask]
        edit_pass = group[read_mask & umi_mask & edit_mask]
        prime_pass = group[read_mask & umi_mask & edit_mask & prime_mask]
        mismatch_pass = group[read_mask & umi_mask & edit_mask & prime_mask & mismatch_mask]
        reference_evidence_pass = group[read_mask & umi_mask & edit_mask & prime_mask & mismatch_mask & reference_evidence]
        reference_pass = group[read_mask & umi_mask & edit_mask & prime_mask & mismatch_mask & reference_mask]
        sequence_pass = group[read_mask & umi_mask & edit_mask & prime_mask & mismatch_mask & reference_mask & sequence_mask]
        conflicting_loci = _conflicting_loci(group, locus_order)
        source_duplicate_conflict = bool(SOURCE_DUPLICATE_CONFLICT_COLUMN in group.columns and _as_bool_series(group[SOURCE_DUPLICATE_CONFLICT_COLUMN]).any())
        selected = _select_barcodes(sequence_pass, locus_order)
        total = len(group)
        read_n = len(read_pass)
        umi_n = len(umi_pass)
        edit_n = len(edit_pass)
        prime_n = len(prime_pass)
        mismatch_n = len(mismatch_pass)
        reference_evidence_n = len(reference_evidence_pass)
        reference_n = len(reference_pass)
        sequence_n = len(sequence_pass)
        used_n = len(selected)
        complexity = used_n
        status = _status(total, read_n, umi_n, edit_n, prime_n, mismatch_n, reference_evidence_n, reference_n, sequence_n, used_n, complexity, qc)
        if source_duplicate_conflict:
            status = 'SOURCE_DUPLICATE_CONFLICT'
        elif conflicting_loci:
            status = 'CONFLICT'
        passed = status == 'PASS'
        loci_text, signature, summary = _summarize_selected(selected, locus_order)
        clone_id = _clone_id(signature) if passed else ''
        min_read_count = int(selected['DARLIN_Read_Count'].min()) if used_n else 0
        min_umi_count = int(selected['DARLIN_UMI_Count'].min()) if used_n else 0
        min_cell_read_fraction = float(selected['DARLIN_Cell_Read_Fraction'].min()) if used_n else 0.0
        read_support_fraction = _support_fraction(selected, sequence_pass, 'read')
        umi_support_fraction = _support_fraction(selected, sequence_pass, 'umi_count')
        reference_qc_status = _reference_qc_status(mismatch_n, reference_evidence_n, reference_n, bool(qc['require_reference_evidence']))
        max_reference_sample_count_observed = _max_reference_sample_count(mismatch_pass, OBSERVED_EVIDENCE_COLUMN)
        max_qc_passing_biological_sample_count = _max_reference_sample_count(mismatch_pass, QC_PASSING_EVIDENCE_COLUMN)
        prime_evidence, max_3prime_mismatches = _summarize_3prime_evidence(selected)
        max_ambiguous_bases = _ambiguous_base_count(selected['DARLIN_Barcode_Value']) if used_n else 0
        double_reasons = []
        if status == 'TOO_MANY_BARCODES':
            double_reasons.append('too_many_darlin_barcodes')
        if used_n and sequence_n > used_n:
            double_reasons.append('multiple_passing_darlin_candidates')
        if conflicting_loci:
            double_reasons.append('conflicting_locus_barcodes:' + ','.join(conflicting_loci))
        aggregation_warnings = []
        duplicate_count = int(_numeric_series(group, SOURCE_DUPLICATE_COUNT_COLUMN).sum())
        if duplicate_count:
            aggregation_warnings.append(f'SOURCE_EXACT_DUPLICATE_DEDUPLICATED:{duplicate_count}')
        if source_duplicate_conflict:
            aggregation_warnings.append('SOURCE_DUPLICATE_CONFLICT')
        rows.append({'Sample_ID': standard_sample_id, 'DARLIN_Clone_ID': clone_id, 'DARLIN_QC_Pass': passed, 'DARLIN_QC_Status': status, 'DARLIN_Loci': loci_text if passed else '', 'DARLIN_Barcode_Summary': summary if passed else '', 'DARLIN_Total_Read_Count': int(selected['DARLIN_Read_Count'].sum()) if used_n else 0, 'DARLIN_Total_UMI_Count': int(selected['DARLIN_UMI_Count'].sum()) if used_n else 0, 'DARLIN_Min_Barcode_Read_Count': min_read_count if passed else 0, 'DARLIN_Min_Barcode_UMI_Count': min_umi_count if passed else 0, 'DARLIN_Min_Cell_Read_Fraction': min_cell_read_fraction if passed else 0.0, 'DARLIN_3Prime_Evidence': prime_evidence if passed else '', 'DARLIN_Max_3Prime_Mismatch_Count': max_3prime_mismatches if passed else 0, 'DARLIN_Max_Ambiguous_Base_Count': max_ambiguous_bases if passed else 0, 'DARLIN_Barcode_Complexity': complexity, 'DARLIN_Assignment_Method': ASSIGNMENT_METHOD, 'DARLIN_Assignment_Read_Support_Fraction': read_support_fraction, 'DARLIN_Assignment_UMI_Support_Fraction': umi_support_fraction, 'DARLIN_Reference_QC_Status': reference_qc_status, 'DARLIN_Reference_Evidence_Method': REFERENCE_EVIDENCE_METHOD, 'DARLIN_Reference_Evidence_Complete': bool(mismatch_n) and reference_evidence_n == mismatch_n, 'DARLIN_Reference_Evidence_Barcodes': reference_evidence_n, 'DARLIN_Max_Reference_Sample_Count_Observed': max_reference_sample_count_observed, 'DARLIN_Max_QC_Passing_Biological_Sample_Count': max_qc_passing_biological_sample_count, 'DARLIN_Double_Cell_Suspect': bool(double_reasons), 'DARLIN_Double_Cell_Reason': ';'.join(double_reasons), 'DARLIN_Aggregation_Warnings': ';'.join(aggregation_warnings), 'DARLIN_Source_Duplicate_Count': duplicate_count, 'DARLIN_Barcodes_Total': total, 'DARLIN_Barcodes_After_Read_QC': read_n, 'DARLIN_Barcodes_After_UMI_QC': umi_n, 'DARLIN_Barcodes_After_Edit_QC': edit_n, 'DARLIN_Barcodes_After_3Prime_QC': prime_n, 'DARLIN_Barcodes_After_3Prime_Mismatch_QC': mismatch_n, 'DARLIN_Barcodes_After_Reference_QC': reference_n, 'DARLIN_Barcodes_After_Sequence_QC': sequence_n, 'DARLIN_Barcodes_Used': used_n if passed else 0, 'QC_Extraction_Read_Cutoff': extraction_read_cutoff, 'QC_Extraction_3Prime_Max_Mismatches': extraction_3prime_max_mismatches, 'QC_Extraction_Min_Cell_Read_Fraction': extraction_min_cell_read_fraction, 'QC_Min_Read_Count': qc['min_read_count'], 'QC_Min_UMI_Count': qc['min_umi_count'], 'QC_Require_Edited_Barcode': qc['require_edited_barcode'], 'QC_Require_3Prime_Detected': qc['require_3prime_detected'], 'QC_Max_3Prime_Mismatches': qc['max_3prime_mismatches'], 'QC_Allow_Read_End_Rescue': qc['allow_read_end_rescue'], 'QC_Min_Read_End_Read_Count': qc['min_read_end_read_count'], 'QC_Min_Read_End_UMI_Count': qc['min_read_end_umi_count'], 'QC_Min_Read_End_Cell_Read_Fraction': qc['min_read_end_cell_read_fraction'], 'QC_Min_Read_End_Barcode_Length': qc['min_read_end_barcode_length'], 'QC_Max_Ambiguous_Bases': qc['max_ambiguous_bases'], 'QC_Max_Barcode_Complexity': qc['max_barcode_complexity'], 'QC_Max_Reference_Sample_Count': qc['max_reference_sample_count'], 'QC_Require_Reference_Evidence': qc['require_reference_evidence'], 'Barcode_Map_SHA256': barcode_map_sha256, 'Barcode_Map_Rows': int(barcode_map_rows)})
    final = pd.DataFrame(rows, columns=FINAL_COLUMNS)
    if not final.empty:
        final = _collapse_duplicate_sample_ids(final)
        final = final.sort_values('Sample_ID').reset_index(drop=True)
    _ensure_parent(output)
    final.to_csv(output, index=False)


def barcode_table_main(argv=None):
    parser = argparse.ArgumentParser(description='Build DARLIN barcode and clone ID tables.')
    subparsers = parser.add_subparsers(dest='command', required=True)
    merge_parser = subparsers.add_parser('merge-barcodes', help='Merge sample barcode calls for one template.')
    merge_parser.add_argument('--input-files', nargs='+', required=True, help='Per-sample DARLIN_Barcodes.csv files.')
    merge_parser.add_argument('--samples', nargs='+', required=True, help='Sample names matching --input-files order.')
    merge_parser.add_argument('--template', required=True, help='Template name represented by the input files.')
    merge_parser.add_argument('--locus', required=True, help='Locus label written to the merged barcode table.')
    merge_parser.add_argument('--configfile', default=None, help='RNA_DARLIN YAML used for non-recurrence QC thresholds.')
    merge_parser.add_argument('--output', required=True, help='Merged barcode CSV path.')
    merge_parser.add_argument('--checksum-output', default=None, help='Optional SHA-256 sidecar for the persistent merged table.')
    clone_parser = subparsers.add_parser('build-clone-table', help='Build final DARLIN_Clone_ID.csv.')
    clone_parser.add_argument('--barcode-files', nargs='+', required=True, help='Merged barcode CSV files.')
    clone_parser.add_argument('--templates', nargs='+', required=True, help='Template names matching --barcode-files order.')
    clone_parser.add_argument('--loci', nargs='+', required=True, help='Locus labels matching --barcode-files order.')
    clone_parser.add_argument('--checksum-files', nargs='+', default=None, help='Optional SHA-256 sidecars matching --barcode-files order.')
    clone_parser.add_argument('--configfile', required=True, help='RNA_DARLIN YAML config.')
    clone_parser.add_argument('--output', required=True, help='Final DARLIN clone table path.')
    args = parser.parse_args(argv)
    if args.command == 'merge-barcodes':
        merge_barcodes(args.input_files, args.samples, args.template, args.locus, args.output, configfile=args.configfile, checksum_output=args.checksum_output)
    elif args.command == 'build-clone-table':
        if len(args.barcode_files) != len(args.templates) or len(args.templates) != len(args.loci):
            raise ValueError('--barcode-files, --templates, and --loci must have the same length')
        build_clone_table(args.barcode_files, args.templates, args.loci, args.configfile, args.output, checksum_files=args.checksum_files)


EXTRACT_OUTPUT_COLUMNS = ['library', 'cell_id', 'cell_bc', 'clone_id', '3prime_detected', '3prime_match_type', '3prime_mismatch_count', 'read', 'umi_count', 'bc_length', 'total_read_per_cell', 'total_read_per_clone', 'cell_read_fraction', 'clone_read_fraction', 'edited']


READ_FATE_SIDECAR_COLUMNS = ['read_id', 'template', 'fate', 'cell_barcode', 'clone_id', 'detail']


RAW_BARCODE_RE = re.compile('[ACGT]+')


def _write_empty_output(output, pandas_module):
    """写出只有表头的固定 clone 表，使空样品仍保留明确 schema。"""
    pandas_module.DataFrame(columns=EXTRACT_OUTPUT_COLUMNS).to_csv(output, index=False)


def _write_read_fate_sidecar(output, template, raw_reads, total_reads, candidate_reads, final_calls, clone_lineage=None):
    """Write one deterministic template-level fate row per raw paired read.

    ``raw_reads`` is either the streamed read-ID inventory path written by the
    multi-template driver or an in-memory inventory DataFrame (unit tests). Candidates join the inventory by global read
    ordinal, so neither the raw read-ID inventory nor read-level evidence must
    stay fully resident, and the FASTQ order is preserved without sorting
    millions of read IDs. Denoise lineage keys include the 3' evidence
    dimensions used by the call table, preventing collisions between otherwise
    identical cell/UMI/input-sequence candidates with different evidence states.
    Fate columns are computed vectorized over candidates and scattered into
    full-length arrays; rows stream to disk in blocks instead of one Python
    ``writerow`` call per read.
    """
    evidence_cols = ['3prime_detected', '3prime_match_type', '3prime_mismatch_count']
    required_candidates = {'read_ordinal', 'cell_bc', 'cell_id', 'umi', 'clone_id', *evidence_cols}
    missing_candidates = required_candidates.difference(candidate_reads.columns) if not candidate_reads.empty else set()
    if missing_candidates:
        raise ValueError(f'DARLIN candidate reads lack read-fate columns: {sorted(missing_candidates)}')
    if candidate_reads.empty:
        candidate_reads = candidate_reads.reindex(columns=sorted(required_candidates))
    total_reads = int(total_reads)
    ordinals = candidate_reads['read_ordinal'].to_numpy(dtype=np.int64)
    if ordinals.size:
        if int(ordinals.min()) < 0 or int(ordinals.max()) >= total_reads:
            raise ValueError('DARLIN candidate read ordinals fall outside the raw read inventory')
        if np.any(np.diff(ordinals) <= 0):
            raise ValueError('DARLIN template parser emitted duplicate candidate read ordinals')
    candidate_cell_barcodes = candidate_reads['cell_bc'].to_numpy(dtype=object, copy=False)
    candidate_input_clones = candidate_reads['clone_id'].to_numpy(dtype=object, copy=False)
    candidate_detected = candidate_reads['3prime_detected'].astype(bool).to_numpy()
    candidate_match_types = candidate_reads['3prime_match_type'].to_numpy(dtype=object, copy=False)
    candidate_mismatches = pd.to_numeric(candidate_reads['3prime_mismatch_count'], errors='coerce').fillna(-1).astype(int).to_numpy()

    # Resolve denoise parents with a hash join; only lineage violations fall
    # back to constructing the offending key for the fail-closed message.
    if clone_lineage is None:
        corrected_clones = candidate_input_clones
    else:
        lineage_cols = ['cell_id', 'cell_bc', 'umi', *evidence_cols, 'input_clone_id', 'corrected_clone_id']
        missing_lineage = set(lineage_cols).difference(clone_lineage.columns)
        if missing_lineage:
            raise ValueError(f'DARLIN denoise lineage lacks columns: {sorted(missing_lineage)}')
        lineage_frame = pd.DataFrame({
            'cell_id': clone_lineage['cell_id'].astype(str),
            'cell_bc': clone_lineage['cell_bc'].astype(str),
            'umi': clone_lineage['umi'].astype(str),
            '3prime_detected': clone_lineage['3prime_detected'].astype(bool),
            '3prime_match_type': clone_lineage['3prime_match_type'].astype(str),
            '3prime_mismatch_count': pd.to_numeric(clone_lineage['3prime_mismatch_count'], errors='coerce').fillna(-1).astype(int),
            'input_clone_id': clone_lineage['input_clone_id'].astype(str),
        })
        lineage_frame['corrected_clone_id'] = clone_lineage['corrected_clone_id'].map(
            lambda value: '' if value is None or value != value else str(value)
        )
        unique_parents = lineage_frame.drop_duplicates()
        conflict_mask = unique_parents.duplicated(subset=lineage_cols[:-1], keep=False)
        if bool(conflict_mask.any()):
            row = unique_parents.loc[unique_parents.index[np.nonzero(conflict_mask.to_numpy())[0][0]]]
            key = (str(row['cell_id']), str(row['cell_bc']), str(row['umi']), bool(row['3prime_detected']), str(row['3prime_match_type']), int(row['3prime_mismatch_count']), str(row['input_clone_id']))
            raise ValueError(f'DARLIN denoise lineage maps one full-evidence candidate to multiple parents: {key!r}')
        candidate_keys = pd.DataFrame({
            'cell_id': candidate_reads['cell_id'].astype(str),
            'cell_bc': candidate_reads['cell_bc'].astype(str),
            'umi': candidate_reads['umi'].astype(str),
            '3prime_detected': candidate_detected,
            '3prime_match_type': candidate_match_types,
            '3prime_mismatch_count': candidate_mismatches,
            'input_clone_id': candidate_input_clones,
        })
        joined = candidate_keys.merge(unique_parents, on=lineage_cols[:-1], how='left', sort=False)
        if bool(joined['corrected_clone_id'].isna().any()):
            position = int(np.nonzero(joined['corrected_clone_id'].isna().to_numpy())[0][0])
            row = joined.iloc[position]
            key = (str(row['cell_id']), str(row['cell_bc']), str(row['umi']), bool(row['3prime_detected']), str(row['3prime_match_type']), int(row['3prime_mismatch_count']), str(row['input_clone_id']))
            raise ValueError(f'DARLIN candidate is absent from the denoise lineage: {key!r}')
        corrected_clones = joined['corrected_clone_id'].to_numpy(dtype=object, copy=False)

    final_pairs = set()
    if final_calls is not None and (not final_calls.empty):
        final_pairs = {(str(row.cell_bc), str(row.clone_id)) for row in final_calls[['cell_bc', 'clone_id']].itertuples(index=False)}
    if final_pairs:
        pairs_frame = pd.DataFrame(final_pairs, columns=['cell_bc', 'clone_id'])
        pair_hit = pd.DataFrame({'cell_bc': candidate_cell_barcodes, 'clone_id': corrected_clones}).merge(
            pairs_frame.assign(_final=True), on=['cell_bc', 'clone_id'], how='left', sort=False
        )['_final'].notna().to_numpy()
    else:
        pair_hit = np.zeros(ordinals.size, dtype=bool)

    fate_candidates = np.where(pair_hit, 'DARLIN_ASSIGNED', 'DARLIN_FAILED_QC').astype(object)
    detail_candidates = np.empty(ordinals.size, dtype=object)
    if clone_lineage is not None:
        failed_detail = np.where(
            corrected_clones == '',
            'primer_and_whitelist_passed_but_failed_read_cutoff',
            'primer_and_whitelist_passed_but_not_final',
        ).astype(object)
    else:
        failed_detail = np.array(['primer_and_whitelist_passed_but_not_final'] * ordinals.size, dtype=object)
    detail_candidates[:] = failed_detail
    remapped = pair_hit & (candidate_input_clones != corrected_clones)
    for position in np.nonzero(remapped)[0]:
        detail_candidates[position] = (
            f'supports_final_template_call_after_directional_denoising;input_clone_id={candidate_input_clones[position]};parent_clone_id={corrected_clones[position]}'
        )
    detail_candidates[pair_hit & ~remapped] = 'supports_final_template_call'

    fate_all = np.full(total_reads, 'DARLIN_UNASSIGNED', dtype=object)
    cell_all = np.full(total_reads, '', dtype=object)
    clone_all = np.full(total_reads, '', dtype=object)
    detail_all = np.full(total_reads, 'no_valid_template_evidence', dtype=object)
    fate_all[ordinals] = fate_candidates
    cell_all[ordinals] = candidate_cell_barcodes
    clone_all[ordinals] = np.where(pair_hit, corrected_clones, candidate_input_clones)
    detail_all[ordinals] = detail_candidates

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with atomic_output(output) as temporary:
        block_size = 1_000_000
        written = 0
        with gzip.open(temporary, 'wt', encoding='utf-8', newline='', compresslevel=1) as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(READ_FATE_SIDECAR_COLUMNS)
            inventory_handle = gzip.open(raw_reads, 'rt', encoding='utf-8') if isinstance(raw_reads, (str, os.PathLike)) else None
            try:
                if inventory_handle is not None:
                    line_iterator = iter(inventory_handle)
                else:
                    line_iterator = iter(raw_reads['read_id'].astype(str))
                while True:
                    block = list(islice(line_iterator, block_size))
                    if not block:
                        break
                    read_ids = [line.rstrip('\n') for line in block]
                    offset = written
                    written += len(read_ids)
                    if written > total_reads:
                        raise ValueError('DARLIN read-ID inventory is longer than the parsed read count')
                    end = offset + len(read_ids)
                    writer.writerows(zip(
                        read_ids,
                        (template,) * len(read_ids),
                        fate_all[offset:end],
                        cell_all[offset:end],
                        clone_all[offset:end],
                        detail_all[offset:end],
                        strict=True,
                    ))
            finally:
                if inventory_handle is not None:
                    inventory_handle.close()
        if written != total_reads:
            raise ValueError('DARLIN read-ID inventory does not match the parsed read count')


def _link_fastq(source, target):
    target = Path(target)
    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to(Path(source).resolve())


def _prepare_fastq_dir(fastq_r1, fastq_r2, sample, work_dir):
    fastq_dir = Path(work_dir) / 'input_fastq'
    fastq_dir.mkdir(parents=True, exist_ok=True)
    _link_fastq(fastq_r1, fastq_dir / f'{sample}_R1.fastq.gz')
    _link_fastq(fastq_r2, fastq_dir / f'{sample}_R2.fastq.gz')
    return fastq_dir


def _template_reference(template):
    """Return the unedited reference sequence for a DARLIN template."""
    if template.startswith('cCARLIN'):
        return CA_CARLIN
    if template.startswith('Tigre'):
        return TA_CARLIN
    if template.startswith('Rosa'):
        return RA_CARLIN
    raise ValueError('template must start with one of: cCARLIN, Tigre, Rosa')


def _resolve_read_cutoff(value, df_unique):
    """Resolve numeric or automatic read cutoff."""
    if str(value) == 'NA':
        return estimate_read_cutoff(df_unique)
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


def _load_cell_barcode_whitelist(path, expected_length):
    """Load the first column of the project-generated RNA whitelist."""
    if expected_length not in {8, 10}:
        raise ValueError(f'DARLIN cell barcode length must be 8 or 10, found {expected_length!r}')
    whitelist_path = Path(path)
    barcodes = []
    try:
        with whitelist_path.open('r', encoding='utf-8', newline='') as handle:
            for line_number, row in enumerate(csv.reader(handle, delimiter='\t'), 1):
                if not row or all((not field for field in row)):
                    raise ValueError(f'DARLIN cell barcode whitelist contains a blank line at {line_number}')
                if len(row) not in {1, 2}:
                    raise ValueError(f'DARLIN cell barcode whitelist must have one or two tab-separated columns; found {len(row)} at line {line_number}')
                barcode = row[0]
                if RAW_BARCODE_RE.fullmatch(barcode) is None:
                    raise ValueError(f'DARLIN cell barcode whitelist must contain raw uppercase A/C/G/T; found {barcode!r} at line {line_number}')
                if len(barcode) != expected_length:
                    raise ValueError(f'DARLIN cell barcode whitelist conflicts with the configured length: expected {expected_length}, found {len(barcode)} for {barcode!r} at line {line_number}')
                barcodes.append(barcode)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f'Cannot read DARLIN cell barcode whitelist {whitelist_path}: {exc}') from exc
    if not barcodes:
        raise ValueError(f'DARLIN cell barcode whitelist is empty: {whitelist_path}')
    if len(set(barcodes)) != len(barcodes):
        raise ValueError(f'DARLIN cell barcode whitelist contains duplicate barcodes: {whitelist_path}')
    return tuple(barcodes)


def _combine_streamed_darlin_candidates(frames):
    """Combine chunk-local read counts without retaining raw FASTQ payload."""
    if not frames:
        return pd.DataFrame(columns=['cell_bc', 'library', 'cell_id', 'umi', 'umi_id', 'clone_id', '3prime_detected', '3prime_match_type', '3prime_mismatch_count', 'read'])
    merged = pd.concat(frames, ignore_index=True)
    if merged.empty:
        return merged
    keys = [column for column in merged.columns if column != 'read']
    return merged.groupby(keys, sort=False, observed=True, dropna=False, as_index=False)['read'].sum()


def _apply_umi_consensus(filtered_reads):
    """Resolve sequencing-error variants within one biological cell+UMI.

    Directional edit-distance-1 clustering occurs before the allele read cutoff,
    so multiple noisy reads from one molecule can jointly support the same
    molecule rather than being discarded as separate low-count sequences.
    After remapping, each ``(cell, UMI, consensus allele)`` molecule is
    represented by one row with summed read support and the deterministic
    best 3-prime evidence signature of its member reads; distinct alleles
    captured by the same UMI stay separate fail-closed candidates.

    Returns ``(consensus_frame, remap)``. ``remap`` maps
    ``(cell_id, cell_bc, umi, observed_clone_id)`` to
    ``(consensus_clone_id, 3prime_detected, 3prime_match_type,
    3prime_mismatch_count)``. Read-level evidence must be remapped with the
    same dictionary before the read-fate sidecar is built, otherwise sidecar
    lineage validation fails closed because candidate input alleles no longer
    match the denoise lineage.
    """
    evidence_columns = ('3prime_detected', '3prime_match_type', '3prime_mismatch_count')
    if filtered_reads is None or filtered_reads.empty:
        return filtered_reads, {}
    # 输入帧来自 chunk 聚合管道，为全新对象；就地变更可避免一整份拷贝。
    df = filtered_reads
    local_keys = [key for key in ('library', 'cell_id', 'cell_bc', 'umi') if key in df.columns]
    observed_clone_ids = df['clone_id'].astype(str)
    corrected = observed_clone_ids.copy()
    remap = {}
    if local_keys:
        group_sizes = df.groupby(local_keys, sort=False, observed=True)['clone_id'].transform('size')
        multi_rows = group_sizes.gt(1)
        if bool(multi_rows.any()):
            for _, indices in df.loc[multi_rows].groupby(local_keys, sort=False, observed=True).groups.items():
                indices = list(indices)
                original_ids = df.loc[indices, 'clone_id'].astype(str)
                # umi_tools directional clustering requires equal-length sequences;
                # the production denoise path partitions by barcode length upstream,
                # so consensus must do the same inside each cell+UMI molecule group.
                consensus_ids = list(original_ids)
                lengths = original_ids.str.len()
                for target_length in sorted(lengths.unique()):
                    positions = [position for position, length in enumerate(lengths) if length == target_length]
                    if len(positions) < 2:
                        continue
                    subset = original_ids.iloc[positions]
                    corrected_ids = [str(value) for value in denoise_sequence(subset, df.loc[subset.index, 'read'])]
                    for position, consensus in zip(positions, corrected_ids, strict=True):
                        consensus_ids[position] = consensus
                corrected.loc[indices] = consensus_ids
                cell_id = str(df.loc[indices[0], 'cell_id']) if 'cell_id' in df.columns else ''
                cell_bc = str(df.loc[indices[0], 'cell_bc']) if 'cell_bc' in df.columns else ''
                umi = str(df.loc[indices[0], 'umi']) if 'umi' in df.columns else ''
                for original, consensus in zip(original_ids, consensus_ids, strict=True):
                    if original != consensus:
                        remap[(cell_id, cell_bc, umi, original)] = consensus
    df['clone_id'] = corrected
    df['_observed_clone_id'] = observed_clone_ids
    molecule_keys = [key for key in ('library', 'cell_id', 'cell_bc', 'umi', 'umi_id', 'clone_id') if key in df.columns]
    has_evidence = all(column in df.columns for column in evidence_columns)
    # Vectorized per-molecule aggregation: one global priority sort picks each
    # molecule's winning evidence row; summed reads and the best 5' edit
    # distance attach via transforms instead of a per-molecule Python loop.
    priority_columns = [column for column in ('3prime_detected', '3prime_mismatch_count', 'read') if column in df.columns]
    priority_ascending = [False, True, False][:len(priority_columns)]
    df = df.sort_values(priority_columns, ascending=priority_ascending, kind='mergesort')
    molecule_groups = [df[column] for column in molecule_keys]
    read_values = pd.to_numeric(df['read'], errors='coerce').fillna(0)
    read_sums = read_values.groupby(molecule_groups).transform('sum')
    five_best = None
    if '5prime_edit_distance' in df.columns:
        five_prime = pd.to_numeric(df['5prime_edit_distance'], errors='coerce')
        five_best = five_prime.where(five_prime.ge(0)).groupby(molecule_groups).transform('min').fillna(-1.0)
    remap_evidence = {}
    if has_evidence:
        winner_detected = df.groupby(molecule_groups, sort=False, observed=True)['3prime_detected'].transform('first').astype(bool).to_numpy()
        winner_match = df.groupby(molecule_groups, sort=False, observed=True)['3prime_match_type'].transform('first').astype(str).to_numpy()
        winner_mismatch = pd.to_numeric(df.groupby(molecule_groups, sort=False, observed=True)['3prime_mismatch_count'].transform('first'), errors='coerce').fillna(-1).astype(int).to_numpy()
        row_detected = df['3prime_detected'].astype(bool).to_numpy()
        row_match = df['3prime_match_type'].astype(str).to_numpy()
        row_mismatch = pd.to_numeric(df['3prime_mismatch_count'], errors='coerce').fillna(-1).astype(int).to_numpy()
        evidence_differs = (
            (row_detected != winner_detected)
            | (row_match != winner_match)
            | (row_mismatch != winner_mismatch)
        )
        cell_ids = df['cell_id'].astype(str).to_numpy() if 'cell_id' in df.columns else np.array([''] * len(df), dtype=object)
        cell_bcs = df['cell_bc'].astype(str).to_numpy() if 'cell_bc' in df.columns else np.array([''] * len(df), dtype=object)
        umis = df['umi'].astype(str).to_numpy() if 'umi' in df.columns else np.array([''] * len(df), dtype=object)
        observed_values = df['_observed_clone_id'].astype(str).to_numpy()
        for position in np.nonzero(evidence_differs)[0]:
            remap_evidence[(cell_ids[position], cell_bcs[position], umis[position], observed_values[position])] = (
                bool(winner_detected[position]), str(winner_match[position]), int(winner_mismatch[position]),
            )
    winners = df.drop_duplicates(subset=molecule_keys, keep='first').copy()
    winners['read'] = read_sums.loc[winners.index].astype(int).to_numpy()
    if five_best is not None:
        winners['5prime_edit_distance'] = five_best.loc[winners.index].astype(int).to_numpy()
    winners = winners.drop(columns=['_observed_clone_id'])
    consensus_frame = winners.reset_index(drop=True)
    final_remap = {}
    for key, consensus_clone in remap.items():
        evidence = remap_evidence.get(key)
        final_remap[key] = (consensus_clone, *(evidence if evidence is not None else (None, None, None)))
    for key, evidence in remap_evidence.items():
        if key not in final_remap:
            final_remap[key] = (key[3], *evidence)
    return consensus_frame, final_remap


def _remap_read_level_alleles(candidates, remap):
    """Apply UMI-consensus clone and evidence remapping to read-level rows.

    ``remap`` is tiny (only multi-allele molecules), so matching uses a hash
    join instead of materializing one Python tuple per candidate read.
    """
    if candidates is None or candidates.empty or not remap:
        return candidates
    required = {'cell_id', 'cell_bc', 'umi', 'clone_id'}
    if not required.issubset(candidates.columns):
        raise ValueError('DARLIN read-level candidates lack UMI-consensus remap columns')
    keys = ['cell_id', 'cell_bc', 'umi', 'clone_id']
    remap_frame = pd.DataFrame(list(remap.keys()), columns=keys)
    remap_frame['_consensus_clone'] = [entry[0] for entry in remap.values()]
    remap_frame['_detected'] = [entry[1] for entry in remap.values()]
    remap_frame['_match_type'] = [entry[2] for entry in remap.values()]
    remap_frame['_mismatch'] = [entry[3] for entry in remap.values()]
    joined = candidates.merge(remap_frame, on=keys, how='left', sort=False)
    hit = joined['_consensus_clone'].notna().to_numpy()
    if not hit.any():
        return candidates
    positions = np.nonzero(hit)[0]
    clone_values = candidates['clone_id'].to_numpy(dtype=object, copy=True)
    detected_values = candidates['3prime_detected'].to_numpy(copy=True) if '3prime_detected' in candidates.columns else None
    match_values = candidates['3prime_match_type'].to_numpy(dtype=object, copy=True) if '3prime_match_type' in candidates.columns else None
    mismatch_values = pd.to_numeric(candidates['3prime_mismatch_count'], errors='coerce').fillna(-1).to_numpy(copy=True) if '3prime_mismatch_count' in candidates.columns else None
    for position in positions:
        entry = remap.get(tuple(joined.loc[position, key] for key in keys))
        if entry is None:
            continue
        clone_values[position] = entry[0]
        if entry[1] is not None and detected_values is not None:
            detected_values[position] = entry[1]
            match_values[position] = entry[2]
            mismatch_values[position] = entry[3]
    candidates['clone_id'] = clone_values
    if detected_values is not None:
        candidates['3prime_detected'] = detected_values
        candidates['3prime_match_type'] = match_values
        candidates['3prime_mismatch_count'] = mismatch_values
    return candidates


def _finalize_darlin_template_calls(filtered_reads, *, template, read_cutoff_setting, min_reads_per_allele_group, min_cell_read_fraction):
    """Apply unchanged biological thresholds after streaming/template parsing."""
    if filtered_reads.empty:
        return pd.DataFrame(columns=EXTRACT_OUTPUT_COLUMNS), None
    filtered_reads = filtered_reads.copy()
    filtered_reads['bc_length'] = filtered_reads['clone_id'].str.len()
    df_unique = unique_cell_cutoff_stats(filtered_reads)
    read_cutoff = _resolve_read_cutoff(read_cutoff_setting, df_unique)
    denoised_by_length = []
    lineage_by_length = []
    for target_length in sorted(filtered_reads['bc_length'].unique()):
        length_group = filtered_reads[filtered_reads['bc_length'] == target_length]
        high_quality, length_lineage = denoise_clonal_data(
            length_group,
            read_cutoff=read_cutoff,
            group_keys=['library', 'cell_id', 'cell_bc', 'clone_id', 'umi', '3prime_detected', '3prime_match_type', '3prime_mismatch_count'],
            return_lineage=True,
        )
        grouped = group_cells(high_quality, group_keys=['library', 'cell_id', 'cell_bc', 'clone_id', '3prime_detected', '3prime_match_type', '3prime_mismatch_count'])
        grouped['bc_length'] = target_length
        denoised_by_length.append(grouped)
        lineage_by_length.append(length_lineage)
    lineage = pd.concat(lineage_by_length, ignore_index=True) if lineage_by_length else None
    if not denoised_by_length:
        return pd.DataFrame(columns=EXTRACT_OUTPUT_COLUMNS), lineage
    denoised = pd.concat(denoised_by_length, ignore_index=True)
    if denoised.empty:
        return pd.DataFrame(columns=EXTRACT_OUTPUT_COLUMNS), lineage
    clone_max_reads = denoised.groupby('clone_id', sort=False)['read'].transform('max')
    clone_cutoff = np.maximum(float(min_reads_per_allele_group), 0.01 * clone_max_reads)
    denoised = denoised.loc[denoised['read'] >= clone_cutoff].copy()
    if denoised.empty:
        return pd.DataFrame(columns=EXTRACT_OUTPUT_COLUMNS), lineage
    denoised = add_read_fractions(denoised)
    final = denoised[denoised['cell_read_fraction'] >= float(min_cell_read_fraction)].copy()
    if final.empty:
        return pd.DataFrame(columns=EXTRACT_OUTPUT_COLUMNS), lineage
    final = check_editing(final, template)
    reference_sequence = _template_reference(template)
    final = final[final['clone_id'] != reference_sequence].copy()
    return final, lineage


def extract_barcodes_multi(args):
    """Parse each FASTQ pair once and evaluate all configured DARLIN templates.

    Bounded-memory streaming contract: read IDs stream to an inventory file,
    read-level candidates stream to per-template files, and only molecule-level
    aggregates stay resident, so peak RSS scales with molecular diversity
    instead of total read count. Read ordinals (FASTQ order) replace read-ID
    sets for the barcode metrics and sidecar row lookup.
    """
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    parent = project_temporary_parent(Path(args.work_dir) / 'stream')
    parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.darlin-stream.', dir=parent) as temporary:
        work_dir = Path(temporary)
        fastq_dir = _prepare_fastq_dir(args.fastq_r1, args.fastq_r2, args.sample, work_dir)
        ref_barcodes = _load_cell_barcode_whitelist(args.cell_barcode_whitelist, args.cell_barcode_length)
        barcode_index = _build_unique_hamming1_barcode_index(ref_barcodes)
        templates = list(dict.fromkeys(args.template))
        if not templates:
            raise ValueError('At least one DARLIN template is required')

        inventory_path = work_dir / 'read_ids.txt.gz'
        candidate_streams = {}
        candidate_writers = {}
        candidate_row_counts = {template: 0 for template in templates}
        filtered_chunks = {template: [] for template in templates}
        use_short_3prime_by_template = {}
        correction_counts = {'exact': 0, 'corrected_hamming1': 0, 'ambiguous': 0, 'unmatched': 0}
        total_reads = 0
        chunk_count = 0
        started = time.time()
        candidate_columns = ['read_ordinal', 'cell_bc', 'cell_id', 'umi', 'clone_id', '3prime_detected', '3prime_match_type', '3prime_mismatch_count']

        with gzip.open(inventory_path, 'wt', encoding='utf-8', newline='', compresslevel=1) as inventory_handle, \
             ExitStack() as candidate_stack:
            for template in templates:
                candidate_path = work_dir / f'candidates_{template}.tsv.gz'
                handle = candidate_stack.enter_context(gzip.open(candidate_path, 'wt', encoding='utf-8', newline='', compresslevel=1))
                candidate_streams[template] = candidate_path
                candidate_writers[template] = csv.writer(handle, delimiter='\t', lineterminator='\n')
                candidate_writers[template].writerow(candidate_columns)

            for raw_chunk in iter_CARLIN_raw_read_chunks(
                str(fastq_dir), args.sample,
                cell_barcode_length=args.cell_barcode_length,
                umi_length=args.umi_length,
                chunk_size=args.chunk_size,
            ):
                chunk_count += 1
                total_reads += len(raw_chunk)
                inventory_handle.write('\n'.join(raw_chunk['read_id'].astype(str).tolist()))
                inventory_handle.write('\n')
                corrected = correct_darlin_cell_barcodes(raw_chunk, barcode_index=barcode_index)
                counts = corrected['cell_bc_correction_status'].value_counts()
                for status in correction_counts:
                    correction_counts[status] += int(counts.get(status, 0))
                valid_chunk = corrected[corrected['cell_bc'].ne('')].drop(columns=['cell_bc_original', 'cell_bc_correction_status'])
                if chunk_count % 20 == 0:
                    rate = total_reads / max(time.time() - started, 1e-9)
                    print(f'[extract-barcodes-multi] sample={args.sample} chunks={chunk_count} reads={total_reads} rate={rate:.0f} reads/s rss={_peak_rss_mb():.0f}MB', flush=True)
                if valid_chunk.empty:
                    continue
                chunk_length = int(valid_chunk['clone_id'].iloc[:100].str.len().mean())
                for template in templates:
                    inferred_short = chunk_length < 300
                    pinned = use_short_3prime_by_template.get(template)
                    if pinned is None:
                        use_short_3prime_by_template[template] = inferred_short
                    elif pinned != inferred_short:
                        raise ValueError(
                            f'DARLIN primer-variant decision changed mid-stream for template {template!r} '
                            f'in sample {args.sample!r}: first chunks implied '
                            f'{"short" if pinned else "long"} 3-prime primer, a later chunk implies '
                            f'{"short" if inferred_short else "long"}; refusing to switch primers on a mixed-length stream'
                        )

                    def collect(message):
                        return None

                    filtered, candidates = CARLIN_preprocessing(
                        valid_chunk,
                        template=template,
                        ref_cell_barcodes=None,
                        seq_3prime_upper_N=args.seq_3prime_upper_n,
                        seq_3prime_max_mismatches=args.seq_3prime_max_mismatches,
                        seq_5prime_max_edits=args.seq_5prime_max_edits,
                        seq_5prime_upper_N=None,
                        return_read_level=True,
                        use_short_3prime=use_short_3prime_by_template[template],
                        progress=collect,
                    )
                    if not filtered.empty:
                        filtered_chunks[template].append(filtered)
                    if not candidates.empty:
                        candidate_writers[template].writerows(
                            zip(
                                candidates['read_ordinal'].to_numpy(),
                                candidates['cell_bc'].astype(str).to_numpy(),
                                candidates['cell_id'].astype(str).to_numpy(),
                                candidates['umi'].astype(str).to_numpy(),
                                candidates['clone_id'].astype(str).to_numpy(),
                                candidates['3prime_detected'].astype(bool).astype(int).to_numpy(),
                                candidates['3prime_match_type'].astype(str).to_numpy(),
                                pd.to_numeric(candidates['3prime_mismatch_count'], errors='coerce').fillna(-1).astype(int).to_numpy(),
                                strict=True,
                            )
                        )
                        candidate_row_counts[template] += len(candidates)

        candidate_ordinal_arrays = {}
        dual_anchor_ordinal_arrays = {}
        for template in templates:
            template_started = time.time()
            output = output_root / template / args.sample / 'DARLIN_Barcodes.csv'
            sidecar = output_root / template / args.sample / 'read_fate.tsv.gz'
            output.parent.mkdir(parents=True, exist_ok=True)
            filtered = _combine_streamed_darlin_candidates(filtered_chunks[template])
            # UMI-level allele consensus runs before the read cutoff so sequencing-error
            # variants of one molecule jointly support that molecule (documented
            # production architecture, docs/contracts/darlin-design.md item 5). Read-level
            # evidence rows are remapped with the same dictionary so sidecar lineage
            # validation keeps matching candidate alleles to the denoise lineage.
            filtered, consensus_remap = _apply_umi_consensus(filtered)
            candidates = _read_candidate_stream(candidate_streams[template], candidate_columns)
            if candidates.empty:
                candidate_ordinal_arrays[template] = np.empty(0, dtype=np.int64)
                dual_anchor_ordinal_arrays[template] = np.empty(0, dtype=np.int64)
            else:
                candidate_ordinal_arrays[template] = candidates['read_ordinal'].to_numpy(dtype=np.int64)
                detected = candidates['3prime_detected'].astype(bool)
                dual_anchor_ordinal_arrays[template] = candidates.loc[detected, 'read_ordinal'].to_numpy(dtype=np.int64)
            candidates = _remap_read_level_alleles(candidates, consensus_remap)
            final, lineage = _finalize_darlin_template_calls(
                filtered,
                template=template,
                read_cutoff_setting=args.read_cutoff,
                min_reads_per_allele_group=args.min_reads_per_allele_group,
                min_cell_read_fraction=args.min_cell_read_fraction,
            )
            if final.empty:
                _write_empty_output(output, pd)
            else:
                final.to_csv(output, index=False)
            _write_read_fate_sidecar(sidecar, template, inventory_path, total_reads, candidates, final, clone_lineage=lineage)
            print(
                f'[extract-barcodes-multi] template={template} sample={args.sample} '
                f'candidate_reads={candidate_row_counts[template]} final_calls={0 if final.empty else len(final)} '
                f'elapsed={time.time() - template_started:.1f}s',
                flush=True,
            )
            # Release per-template payload before the next template's readback.
            filtered_chunks[template] = []
            del candidates, filtered, final, lineage, consensus_remap

        candidate_ordinal_union = np.unique(np.concatenate(list(candidate_ordinal_arrays.values()))) if candidate_ordinal_arrays else np.empty(0, dtype=np.int64)
        dual_anchor_ordinal_union = np.unique(np.concatenate(list(dual_anchor_ordinal_arrays.values()))) if dual_anchor_ordinal_arrays else np.empty(0, dtype=np.int64)
        metrics_path = Path(args.barcode_metrics)
        metrics_path.parent.mkdir(parents=True, exist_ok=True)
        accepted = correction_counts['exact'] + correction_counts['corrected_hamming1']
        with metrics_path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['sample', 'source_type', 'total_reads', 'exact_barcode_reads', 'corrected_hamming1_reads', 'ambiguous_barcode_reads', 'unmatched_barcode_reads', 'accepted_barcode_reads', 'accepted_fraction', 'correction_fraction_of_accepted', 'darlin_candidate_reads_any_template', 'dual_anchor_reads_any_template'])
            writer.writerow([
                args.sample, 'rna_rescue' if str(args.sample).endswith('__RNA') else 'targeted',
                total_reads, correction_counts['exact'], correction_counts['corrected_hamming1'],
                correction_counts['ambiguous'], correction_counts['unmatched'], accepted,
                f'{accepted / total_reads if total_reads else 0.0:.8f}',
                f'{correction_counts["corrected_hamming1"] / accepted if accepted else 0.0:.8f}',
                int(candidate_ordinal_union.size), int(dual_anchor_ordinal_union.size),
            ])



def _read_candidate_stream(path, columns):
    """Read back one streamed per-template candidate table."""
    frame = pd.read_csv(
        path, sep='\t', compression='gzip', dtype={'cell_bc': str, 'cell_id': str, 'umi': str, 'clone_id': str, '3prime_match_type': str},
    )
    if frame.empty:
        return pd.DataFrame(columns=columns)
    frame['3prime_detected'] = frame['3prime_detected'].astype(bool)
    return frame


def _peak_rss_mb():
    """Return peak resident set size in MB for progress reporting."""
    try:
        import resource
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # macOS reports ru_maxrss in bytes, Linux in kilobytes.
        return value / (1024.0 * 1024.0) if sys.platform == 'darwin' else value / 1024.0
    except Exception:
        return float('nan')


def extract_barcodes_multi_main(argv=None):
    parser = argparse.ArgumentParser(description='Single-pass multi-template DARLIN extraction.')
    parser.add_argument('--sample', required=True)
    parser.add_argument('--template', action='append', required=True, help='Repeat for each DARLIN template.')
    parser.add_argument('--fastq-r1', required=True)
    parser.add_argument('--fastq-r2', required=True)
    parser.add_argument('--cell-barcode-whitelist', required=True)
    parser.add_argument('--cell-barcode-length', type=int, choices=(8, 10), required=True)
    parser.add_argument('--umi-length', type=int, choices=(8,), default=8)
    parser.add_argument('--output-root', required=True, help='DARLIN work root containing template/sample outputs.')
    parser.add_argument('--barcode-metrics', required=True, help='Per-sample barcode correction metrics TSV.')
    parser.add_argument('--work-dir', required=True)
    parser.add_argument('--chunk-size', type=int, default=100000)
    parser.add_argument('--min-reads-per-allele-group', type=int, required=True)
    parser.add_argument('--seq-5prime-max-edits', type=int, default=1)
    parser.add_argument('--seq-3prime-upper-n', type=int, required=True)
    parser.add_argument('--seq-3prime-max-mismatches', type=int, default=0)
    parser.add_argument('--min-cell-read-fraction', type=float, default=0.8)
    parser.add_argument('--read-cutoff', default='NA')
    extract_barcodes_multi(parser.parse_args(argv))
