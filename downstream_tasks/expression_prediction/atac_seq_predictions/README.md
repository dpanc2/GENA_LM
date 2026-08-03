# ATAC-seq predictions

This workflow predicts ATAC-seq signal for 10 kb regions from BED files using
HEK293 or HCT116 descriptions. The reported score is the mean prediction across
tokens overlapping the central 50 bp (center ±25 bp).

`atac_seq_workflow.ipynb` is a example based on Egor notebook. Batch predictions are run with `predict_atac_center50.py` through `run_atac_center50.sh`.

## Run

Activate the API environment:

```bash
conda activate api
```

Test on 10 intervals:

```bash
GENA_CONDITION=HCT116 GENA_LIMIT=10 \
bash downstream_tasks/expression_prediction/atac_seq_predictions/run_atac_center50.sh
```

Run all intervals in the background:

```bash
GENA_CONDITION=HCT116 \
GENA_OUTPUT_DIR=/home/jovyan/dpanc/komarek/atac_center50_predictions_HCT116 \
nohup bash downstream_tasks/expression_prediction/atac_seq_predictions/run_atac_center50.sh \
> atac_center50_hct116.log 2>&1 &
```

Use `GENA_CONDITION=HEK293` for HEK293 predictions. Output files are
tab-separated and contain one central 50 bp score per BED interval.