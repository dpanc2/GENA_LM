#!/usr/bin/env python3
"""Predict token-level ATAC signal for sliding windows across a genome."""

import argparse
import os
import sys
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--bigwig", required=True, type=Path)
    parser.add_argument("--chrom-sizes", required=True, type=Path)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--asset-root", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--condition", default="HEK293")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--sequence-size", type=int, default=10_000)
    parser.add_argument("--center-bp", type=int, default=1_000)
    parser.add_argument("--chunk-size", type=int, default=20_000)
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def read_bed(path, limit=None):
    records = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 BED columns")
            start, end = int(fields[1]), int(fields[2])
            records.append(
                {
                    "id": fields[3] if len(fields) > 3 and fields[3] else f"{fields[0]}:{start}-{end}",
                    "chrom": fields[0],
                    "bed_start": start,
                    "bed_end": end,
                    "strand": fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else "+",
                }
            )
            if limit is not None and len(records) >= limit:
                break
    return records


def read_chrom_sizes(path):
    sizes = {}
    with path.open() as handle:
        for line in handle:
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
    return sizes


def main():
    args = arguments()
    repo = args.repo_root.expanduser().resolve()
    assets = args.asset_root.expanduser().resolve()
    task_dir = repo / "downstream_tasks" / "expression_prediction"
    api_src = task_dir / "api" / "src"
    sys.path.insert(0, str(api_src))

    os.environ["GENA_REPO_ROOT"] = str(repo)
    os.environ["GENA_DECODER"] = str(assets / "models" / "decoder")

    from gena_expression.conditions import DescriptionLookup
    from gena_expression.contexts import Genome
    from gena_expression.inference import SequenceModel
    import pyBigWig

    required = [
        args.bed,
        args.genome,
        args.chrom_sizes,
        assets / "models" / "full_model" / "pytorch_model.bin",
        assets / "models" / "decoder",
        repo / "data" / "tokenizers" / "t2t_1000h_multi_32k",
        task_dir / "atac_seq_predictions" / "inference_portable.yaml",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required paths:\n" + "\n".join(missing))

    records = read_bed(args.bed, args.limit)
    print(f"Loaded {len(records)} intervals from {args.bed}", flush=True)

    genome = Genome(fasta_path=args.genome)
    model = SequenceModel.load(
        model_cls=f"{task_dir / 'expression_model_final.py'}::ExpressionCounts",
        checkpoint=assets / "models" / "full_model" / "pytorch_model.bin",
        config=task_dir / "atac_seq_predictions" / "inference_portable.yaml",
        dna_tokenizer=repo / "data" / "tokenizers" / "t2t_1000h_multi_32k",
        description_tokenizer="Qwen/Qwen3-Embedding-0.6B",
        dna_max_seq_len=1024,
        num_before=512,
        desc_max_seq_len=510,
        token_len_for_fetch=15,
        device=args.device,
    )
    descriptions = DescriptionLookup(
        json_dir=task_dir / "api" / "data" / "atacseq_generated_descriptions",
        keys=[args.condition],
    )
    condition = descriptions[args.condition]
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)

    args.bigwig.parent.mkdir(parents=True, exist_ok=True)

    for chunk_start in range(0, len(records), args.chunk_size):
        chunk_end = min(chunk_start + args.chunk_size, len(records))
        part_path = args.bigwig.with_name(
            f"{args.bigwig.stem}.intervals_{chunk_start + 1:06d}_{chunk_end:06d}.bw"
        )

        if part_path.exists():
            try:
                old_bigwig = pyBigWig.open(str(part_path))
                is_complete = old_bigwig.header()["nBasesCovered"] > 0
                old_bigwig.close()
            except (RuntimeError, TypeError):
                is_complete = False
            if is_complete:
                print(f"Skipping completed {part_path}", flush=True)
                continue

        print(f"Writing {part_path}", flush=True)
        bigwig = pyBigWig.open(str(part_path), "w")
        bigwig.addHeader(list(chrom_sizes.items()))
        last_end_by_chrom = {}

        try:
            for index in range(chunk_start, chunk_end):
                record = records[index]
                genomic_center = (record["bed_start"] + record["bed_end"]) // 2
                sequence = genome.sequence_around(
                    center=genomic_center,
                    chrom=record["chrom"],
                    strand=record["strand"],
                    size=args.sequence_size,
                )
                prediction = model.predict_sequence(
                    sequence=sequence,
                    center=args.sequence_size // 2,
                    condition=condition,
                )
                frame = prediction.track().to_frame()
                window_start = genomic_center - args.sequence_size // 2
                window_end = window_start + args.sequence_size
                center_start = args.sequence_size // 2 - args.center_bp // 2
                center_end = center_start + args.center_bp

                center_tokens = frame[
                    (frame["end"] > center_start) & (frame["start"] < center_end)
                ]
                bigwig_rows = []

                for token in center_tokens.itertuples(index=False):
                    token_start = max(int(token.start), center_start)
                    token_end = min(int(token.end), center_end)
                    if record["strand"] == "+":
                        genomic_start = window_start + token_start
                        genomic_end = window_start + token_end
                    else:
                        genomic_start = window_end - token_end
                        genomic_end = window_end - token_start

                    genomic_start = max(0, genomic_start)
                    genomic_end = min(chrom_sizes[record["chrom"]], genomic_end)
                    if genomic_end > genomic_start:
                        bigwig_rows.append(
                            (genomic_start, genomic_end, float(token.track))
                        )

                clean_rows = []
                for genomic_start, genomic_end, value in sorted(bigwig_rows):
                    genomic_start = max(
                        genomic_start,
                        last_end_by_chrom.get(record["chrom"], 0),
                    )
                    if genomic_end <= genomic_start:
                        continue
                    clean_rows.append((genomic_start, genomic_end, value))
                    last_end_by_chrom[record["chrom"]] = genomic_end

                if clean_rows:
                    bigwig.addEntries(
                        [record["chrom"]] * len(clean_rows),
                        [row[0] for row in clean_rows],
                        ends=[row[1] for row in clean_rows],
                        values=[row[2] for row in clean_rows],
                    )

                completed = index + 1
                if completed == 1 or completed % 100 == 0:
                    print(f"{completed}/{len(records)}", flush=True)
        finally:
            bigwig.close()

        print(f"Saved {part_path}", flush=True)


if __name__ == "__main__":
    main()
