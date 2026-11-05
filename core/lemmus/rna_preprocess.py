"""DD-MET5 RNA：有界分块并行、原生接头匹配和按输入顺序提交。"""

from __future__ import annotations
from .io import project_temporary_parent
import argparse, io, os, time, json, tempfile, multiprocessing
from pathlib import Path
from collections import deque
from concurrent.futures import ProcessPoolExecutor

_ADAPTERS = ()
_DNA = frozenset("ACGT")
_KEYS = (
    "input_pairs",
    "short_cb_umi",
    "invalid_cb",
    "invalid_umi",
    "short_r2",
    "passed_pairs",
    "adapter_trimmed_pairs",
)
CHUNK_BYTES = 4 * 1024 * 1024


def _initialize(adapters):
    global _ADAPTERS
    from cutadapt.adapters import BackAdapter, RightmostFrontAdapter

    _ADAPTERS = tuple(
        BackAdapter(seq, max_errors=0.1, min_overlap=3)
        if end == "3"
        else RightmostFrontAdapter(seq, max_errors=0.1, min_overlap=3)
        for seq, end in adapters
    )


def _process(payload):
    import dnaio
    from isal import igzip

    started = time.perf_counter()
    cpu = time.process_time()
    with dnaio.open(
        io.BytesIO(payload[0]), io.BytesIO(payload[1]), fileformat="fastq"
    ) as reader:
        pairs = list(reader)
    parsed = time.perf_counter()
    metrics = dict.fromkeys(_KEYS, 0)
    hits = [0] * len(_ADAPTERS)
    kept = []
    for first, second in pairs:
        metrics["input_pairs"] += 1
        if len(first) < 29:
            metrics["short_cb_umi"] += 1
            continue
        sequence = first.sequence
        if not _DNA.issuperset(sequence[:17]):
            metrics["invalid_cb"] += 1
            continue
        if not _DNA.issuperset(sequence[17:29]):
            metrics["invalid_umi"] += 1
            continue
        trimmed = False
        for i, adapter in enumerate(_ADAPTERS):
            match = adapter.match_to(second.sequence)
            if match is not None:
                second = match.trimmed(second)
                hits[i] += 1
                trimmed = True
        metrics["adapter_trimmed_pairs"] += int(trimmed)
        if len(second) < 60:
            metrics["short_r2"] += 1
            continue
        kept.append((first[:29], second))
        metrics["passed_pairs"] += 1
    matched = time.perf_counter()
    a, b = io.BytesIO(), io.BytesIO()
    with dnaio.open(a, b, mode="w", fileformat="fastq") as writer:
        for first, second in kept:
            writer.write(first, second)
    serialized = time.perf_counter()
    out = (
        igzip.compress(a.getvalue(), compresslevel=1, mtime=0),
        igzip.compress(b.getvalue(), compresslevel=1, mtime=0),
    )
    finished = time.perf_counter()
    return (
        out,
        metrics,
        hits,
        {
            "parse_wall_seconds": parsed - started,
            "matching_filtering_wall_seconds": matched - parsed,
            "serialization_wall_seconds": serialized - matched,
            "compression_wall_seconds": finished - serialized,
            "worker_CPU_seconds": time.process_time() - cpu,
        },
    )


def preprocess_droplet_rna(
    r1: Path, r2: Path, output_r1: Path, output_r2: Path, qc: Path, *, threads: int
) -> dict:
    """按原生四步匹配处理配对块；worker 全部成功后提交 FASTQ，QC 最后提交。"""
    from .config import DROPLET_R2_ADAPTERS, DROPLET_SOURCE_COMMIT

    adapters = DROPLET_R2_ADAPTERS
    source_commit = DROPLET_SOURCE_COMMIT
    import dnaio, cutadapt, isal
    from isal import igzip, igzip_threaded

    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")
    destinations = (output_r1, output_r2, qc)
    for path in destinations:
        path.parent.mkdir(parents=True, exist_ok=True)
        project_temporary_parent(path)
    started = time.perf_counter()
    metrics = dict.fromkeys(_KEYS, 0)
    hits = [0] * len(adapters)
    phase = {}
    decode_threads = min(2, max(0, threads - 2)) if threads > 1 else 0
    workers = threads - 1 - decode_threads if threads > 1 else 0
    limit = max(1, 2 * workers)
    chunks = 0
    input_time = output_time = 0.0
    with tempfile.TemporaryDirectory(
        prefix=".droplet-trim.", dir=project_temporary_parent(output_r1)
    ) as tmp:
        a, b = Path(tmp) / "R1.fq.gz", Path(tmp) / "R2.fq.gz"
        with (
            (
                igzip_threaded.open(r1, "rb", threads=1)
                if decode_threads >= 1
                else igzip.open(r1, "rb")
            ) as first,
            (
                igzip_threaded.open(r2, "rb", threads=1)
                if decode_threads >= 2
                else igzip.open(r2, "rb")
            ) as second,
            a.open("wb") as out1,
            b.open("wb") as out2,
        ):

            def consume(result):
                nonlocal output_time
                output, counts, adapter_hits, times = result
                stamp = time.perf_counter()
                out1.write(output[0])
                out2.write(output[1])
                output_time += time.perf_counter() - stamp
                for key, value in counts.items():
                    metrics[key] += value
                for i, value in enumerate(adapter_hits):
                    hits[i] += value
                for key, value in times.items():
                    phase[key] = phase.get(key, 0.0) + value

            source = iter(
                dnaio.read_paired_chunks(first, second, buffer_size=CHUNK_BYTES)
            )

            def next_payload():
                nonlocal input_time, chunks
                stamp = time.perf_counter()
                try:
                    p = next(source)
                except StopIteration:
                    input_time += time.perf_counter() - stamp
                    return None
                p = (bytes(p[0]), bytes(p[1]))
                input_time += time.perf_counter() - stamp
                chunks += 1
                return p

            if workers:
                with ProcessPoolExecutor(
                    max_workers=workers,
                    mp_context=multiprocessing.get_context("spawn"),
                    initializer=_initialize,
                    initargs=(adapters,),
                ) as pool:
                    pending = deque()
                    while True:
                        payload = next_payload()
                        if payload is None:
                            break
                        pending.append(pool.submit(_process, payload))
                        if len(pending) >= limit:
                            consume(pending.popleft().result())
                    for future in pending:
                        consume(future.result())
            else:
                _initialize(adapters)
                while True:
                    payload = next_payload()
                    if payload is None:
                        break
                    consume(_process(payload))
            if not chunks:
                out1.write(igzip.compress(b"", compresslevel=1, mtime=0))
                out2.write(igzip.compress(b"", compresslevel=1, mtime=0))
        if metrics["input_pairs"] != sum(
            metrics[k]
            for k in (
                "short_cb_umi",
                "invalid_cb",
                "invalid_umi",
                "short_r2",
                "passed_pairs",
            )
        ):
            raise ValueError("DD-MET5 preprocessing pair accounting failed")
        metrics.update(
            protocol="droplet",
            chemistry="DD-MET5",
            cutadapt_version=cutadapt.__version__,
            source_commit=source_commit,
            adapters=list(adapters),
            adapter_hits=hits,
            max_error_rate=0.1,
            min_overlap=3,
            min_r2_length=60,
            elapsed_seconds=time.perf_counter() - started,
            execution={
                "threads": threads,
                "worker_processes": workers,
                "background_decoder_threads": decode_threads,
                "max_inflight_chunks": limit,
                "chunk_bytes_per_mate": CHUNK_BYTES,
                "chunks": chunks,
                "isal_version": isal.__version__,
                "input_read_decode_chunk_wall_seconds": input_time,
                "ordered_output_write_wall_seconds": output_time,
                "sum_worker_phase_seconds": phase,
                "phase_scope": "Worker wall sums overlap across processes; matching includes filtering. Input phase includes I/O, decoding wait and chunk copies, not isolated decoder CPU.",
            },
        )
        candidate = Path(tmp) / "qc.json"
        candidate.write_text(json.dumps(metrics, indent=2) + "\n")
        qc.unlink(missing_ok=True)
        os.replace(a, output_r1)
        os.replace(b, output_r2)
        os.replace(candidate, qc)
    return metrics


def preprocess_droplet_main(argv=None):
    """接收执行器分配的总 CPU 预算，不在进程内扩张资源。"""
    parser = argparse.ArgumentParser()
    for arg in ("r1", "r2", "output-r1", "output-r2", "qc"):
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args(argv)
    preprocess_droplet_rna(
        args.r1, args.r2, args.output_r1, args.output_r2, args.qc, threads=args.threads
    )
    return 0
