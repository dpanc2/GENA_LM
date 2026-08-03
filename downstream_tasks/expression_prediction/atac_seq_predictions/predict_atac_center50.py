#!/usr/bin/env python3
"""Predict ATAC signal and report the mean for tokens overlapping the central 50 bp."""

import argparse
import csv
import os
import sys
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--asset-root", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--condition", default="HEK293")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--sequence-size", type=int, default=10_000)
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

    required = [
        args.bed,
        args.genome,
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

    args.output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "id", "chrom", "bed_start", "bed_end", "strand", "sequence_center",
        "n_center_tokens", "predicted_start", "predicted_end", "center_50bp_mean",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for index, record in enumerate(records, 1):
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
            center = args.sequence_size // 2
            center_tokens = frame[
                (frame["end"] > center - 25) & (frame["start"] < center + 25)
            ]
            writer.writerow(
                {
                    **record,
                    "sequence_center": genomic_center,
                    "n_center_tokens": len(center_tokens),
                    "predicted_start": int(frame["start"].min()),
                    "predicted_end": int(frame["end"].max()),
                    "center_50bp_mean": float(center_tokens["track"].mean()),
                }
            )
            handle.flush()
            if index == 1 or index % 100 == 0:
                print(f"{index}/{len(records)}", flush=True)

    print(f"Saved {args.output}", flush=True)


if __name__ == "__main__":
    main()
