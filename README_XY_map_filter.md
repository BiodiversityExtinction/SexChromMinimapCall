# Sex-Chromosome Scaffold Mapping and Classification

This script maps query scaffolds to a reference containing two sex chromosomes (for example `X/Y` or `Z/W`) and classifies scaffolds based on non-overlapping aligned fraction.

Scripts:
- `xy_map_filter.sh` (main implementation)
- `SexChromMinimapCall.sh` (neutral wrapper name)

## What it does

`xy_map_filter.sh` supports two modes:

1. `map + filter` (default)
- Runs `minimap2` to produce a PAF file.
- Filters PAF hits by mapping quality and alignment length.
- Computes non-overlapping query coverage on `chr1` and `chr2` targets.
- Outputs scaffold-level fractions and a call.

2. `--filter-only`
- Skips mapping.
- Uses an existing PAF file for classification.

## Why non-overlapping coverage?

PAF can contain overlapping alignments for the same scaffold and target.
Summing raw alignment lengths can overcount and inflate fractions above 1.
This script merges overlapping query intervals per scaffold per target before computing coverage fractions.

## Requirements

- `bash`
- `awk`
- `sort`
- `minimap2` (only required if not using `--filter-only`)

## Usage

### 1) Filter only from existing PAF

```bash
./xy_map_filter.sh \
  --filter-only \
  --paf Asia_black-Brown_xy.paf \
  --chr1 NC_079873.1 \
  --chr2 NC_079874.1 \
  --label1 X \
  --label2 Y \
  --out scaffold_XY_calls.tsv
```

### 2) Map then filter

```bash
./xy_map_filter.sh \
  --query GCA_009660055.1_ASM966005v1_genomic_1Mb.fna \
  --target /path/to/XY_reference.fna \
  --paf Asia_black-Brown_xy.paf \
  --chr1 NC_079873.1 \
  --chr2 NC_079874.1 \
  --label1 X \
  --label2 Y \
  --out scaffold_XY_calls.tsv
```

## Key parameters

Required:
- `--chr1`: target sequence name for chromosome 1 (must match PAF target column 6 exactly)
- `--chr2`: target sequence name for chromosome 2
- `--label1`: label used for chr1 in calls/output column names (default: `X`)
- `--label2`: label used for chr2 in calls/output column names (default: `Y`)

Mapping inputs (unless `--filter-only`):
- `--query`: query FASTA
- `--target`: target FASTA containing X and Y

Filtering thresholds:
- `--min-mapq` (default `30`): keep alignments with MAPQ >= this value
- `--min-aln` (default `55000`): keep alignments with PAF `$11` >= this value
- `--min-frac` (default `0.8`): base fraction threshold for confident calls
- `--chr1-frac-threshold` (default `--min-frac`): explicit chr1 threshold
- `--chr2-frac-threshold` (default `--min-frac`): explicit chr2 threshold
- `--dominance` (default `2.0`): enrichment ratio for mixed chr1+chr2 hits

Minimap:
- `--threads` (default `16`)
- `--preset` (default `asm10`)
- `--secondary yes|no` (default `no`)

## Output format

Output columns:

- `scaffold`
- `scaffold_len`
- `<label1>_cov_bp`
- `<label2>_cov_bp`
- `<label1>_frac`
- `<label2>_frac`
- `<label1>_hits`
- `<label2>_hits`
- `call`

`<label>_hits` are counts of passing PAF alignment records (after filters), not merged blocks.

## Call logic

- `<label1>_only`: chr2 coverage is 0 and `chr1_frac >= chr1-threshold`
- `<label2>_only`: chr1 coverage is 0 and `chr2_frac >= chr2-threshold`
- `<label1>_enriched`: both present, and `chr1_frac >= chr1-threshold` and `chr1_frac >= chr2_frac * dominance`
- `<label2>_enriched`: both present, and `chr2_frac >= chr2-threshold` and `chr2_frac >= chr1_frac * dominance`
- `-`: all other cases (low confidence, mixed, or below threshold)

## Practical notes

- If chr1/chr2 accessions are swapped, calls will invert. Always verify `--chr1`/`--chr2` labels first.
- `--dominance` is especially useful for pseudoautosomal/shared homologous regions to avoid over-calling.
- For strict assignment, use high thresholds (for example `--chr1-frac-threshold 0.8 --chr2-frac-threshold 0.8`).

## Quick sanity checks

Check header + first rows:

```bash
head -n 10 scaffold_XY_calls.tsv
```

Count confident X or Y calls:

```bash
awk 'NR>1 && ($9=="X_only" || $9=="X_enriched") {x++}
     NR>1 && ($9=="Y_only" || $9=="Y_enriched") {y++}
     END{print "X_confident="x+0, "Y_confident="y+0}' scaffold_XY_calls.tsv
```
