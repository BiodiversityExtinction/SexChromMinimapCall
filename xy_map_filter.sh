#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
xy_map_filter.sh: map scaffolds to a sex-chromosome reference (minimap2) and classify chr1/chr2 scaffolds from PAF.

Modes:
  1) Map + filter (default)
  2) Filter only from existing PAF (--filter-only)

Usage:
  xy_map_filter.sh --chr1 TARGET1 --chr2 TARGET2 [options]

Required:
  --chr1 NAME              Target sequence name for chromosome 1 (PAF target column 6)
  --chr2 NAME              Target sequence name for chromosome 2 (PAF target column 6)
  --label1 STR             Label for chr1 in outputs/calls (default: X)
  --label2 STR             Label for chr2 in outputs/calls (default: Y)
  --x NAME                 Backward-compatible alias for --chr1
  --y NAME                 Backward-compatible alias for --chr2

Mapping inputs (required unless --filter-only):
  --query FILE             Query assembly FASTA
  --target FILE            Target reference FASTA containing chr1 and chr2

PAF / outputs:
  --paf FILE               PAF path (output path for mapping; input path for filtering)
                           Default: xy_map.paf
  --out FILE               Filtered scaffold call table
                           Default: scaffold_XY_calls.tsv

Filtering thresholds:
  --min-mapq INT           Keep hits with MAPQ >= this value (default: 30)
  --min-aln INT            Keep hits with alignment block length ($11) >= this value (default: 55000)
  --min-frac FLOAT         Minimum covered fraction to call <label>_only / <label>_enriched (default: 0.8)
  --dominance FLOAT        Ratio for enriched calls when both chr1 and chr2 have hits.
                           <label1>_enriched if chr1_frac >= chr2_frac * dominance (and chr1_frac >= threshold).
                           <label2>_enriched if chr2_frac >= chr1_frac * dominance (and chr2_frac >= threshold). Default: 2.0
  --chr1-frac-threshold FLOAT Explicit threshold to flag scaffold as chr1-like (default: --min-frac)
  --chr2-frac-threshold FLOAT Explicit threshold to flag scaffold as chr2-like (default: --min-frac)
  --x-frac-threshold FLOAT    Backward-compatible alias for --chr1-frac-threshold
  --y-frac-threshold FLOAT    Backward-compatible alias for --chr2-frac-threshold

Minimap options:
  --threads INT            minimap2 threads (default: 16)
  --preset STR             minimap2 preset (default: asm10)
  --secondary yes|no       keep secondary alignments (default: no)

Other:
  --filter-only            Skip minimap2 and only run filtering on --paf
  -h, --help               Show this help

Output columns:
  scaffold  scaffold_len  <label1>_cov_bp  <label2>_cov_bp  <label1>_frac  <label2>_frac  <label1>_hits  <label2>_hits  call

Notes:
  - Coverage is computed as non-overlapping query-span coverage per scaffold per target (chr1/chr2).
  - This avoids double-counting overlapping alignments that can inflate fractions above 1.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }

# --------------------
# Defaults
# --------------------
CHR1_NAME=""
CHR2_NAME=""
LABEL1="X"
LABEL2="Y"
QUERY=""
TARGET=""
PAF="xy_map.paf"
OUT="scaffold_XY_calls.tsv"
THREADS=16
PRESET="asm10"
SECONDARY="no"
FILTER_ONLY=0
MIN_MAPQ=30
MIN_ALN=55000
MIN_FRAC=0.8
DOM=2.0
CHR1_FRAC_THR=""
CHR2_FRAC_THR=""

# --------------------
# Parse args
# --------------------
need_arg() {
  local opt="$1"
  local val="${2-}"
  [[ -n "$val" && "$val" != --* ]] || die "$opt requires a value"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chr1|--x) need_arg "$1" "${2-}"; CHR1_NAME="$2"; shift 2 ;;
    --chr2|--y) need_arg "$1" "${2-}"; CHR2_NAME="$2"; shift 2 ;;
    --label1) need_arg "$1" "${2-}"; LABEL1="$2"; shift 2 ;;
    --label2) need_arg "$1" "${2-}"; LABEL2="$2"; shift 2 ;;
    --query) need_arg "$1" "${2-}"; QUERY="$2"; shift 2 ;;
    --target) need_arg "$1" "${2-}"; TARGET="$2"; shift 2 ;;
    --paf) need_arg "$1" "${2-}"; PAF="$2"; shift 2 ;;
    --out) need_arg "$1" "${2-}"; OUT="$2"; shift 2 ;;
    --threads) need_arg "$1" "${2-}"; THREADS="$2"; shift 2 ;;
    --preset) need_arg "$1" "${2-}"; PRESET="$2"; shift 2 ;;
    --secondary) need_arg "$1" "${2-}"; SECONDARY="$2"; shift 2 ;;
    --filter-only) FILTER_ONLY=1; shift ;;
    --min-mapq) need_arg "$1" "${2-}"; MIN_MAPQ="$2"; shift 2 ;;
    --min-aln) need_arg "$1" "${2-}"; MIN_ALN="$2"; shift 2 ;;
    --min-frac) need_arg "$1" "${2-}"; MIN_FRAC="$2"; shift 2 ;;
    --dominance) need_arg "$1" "${2-}"; DOM="$2"; shift 2 ;;
    --chr1-frac-threshold|--x-frac-threshold) need_arg "$1" "${2-}"; CHR1_FRAC_THR="$2"; shift 2 ;;
    --chr2-frac-threshold|--y-frac-threshold) need_arg "$1" "${2-}"; CHR2_FRAC_THR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

[[ -n "$CHR1_NAME" && -n "$CHR2_NAME" ]] || { usage; die "--chr1 and --chr2 (or --x/--y aliases) are required"; }
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || die "--threads must be a positive integer"
[[ "$MIN_MAPQ" =~ ^[0-9]+$ ]] || die "--min-mapq must be an integer >= 0"
[[ "$MIN_ALN" =~ ^[0-9]+$ ]] || die "--min-aln must be an integer >= 0"
[[ "$SECONDARY" == "yes" || "$SECONDARY" == "no" ]] || die "--secondary must be yes|no"
if [[ -z "$CHR1_FRAC_THR" ]]; then CHR1_FRAC_THR="$MIN_FRAC"; fi
if [[ -z "$CHR2_FRAC_THR" ]]; then CHR2_FRAC_THR="$MIN_FRAC"; fi

if [[ "$FILTER_ONLY" -eq 0 ]]; then
  [[ -n "$QUERY" && -n "$TARGET" ]] || die "--query and --target are required unless --filter-only"
  [[ -r "$QUERY" ]] || die "Query FASTA not readable: $QUERY"
  [[ -r "$TARGET" ]] || die "Target FASTA not readable: $TARGET"
  command -v minimap2 >/dev/null 2>&1 || die "minimap2 not found in PATH"
fi

[[ -r "$PAF" || "$FILTER_ONLY" -eq 0 ]] || die "PAF not readable: $PAF"

command -v awk >/dev/null 2>&1 || die "awk not found in PATH"
command -v sort >/dev/null 2>&1 || die "sort not found in PATH"

echo "[INFO] chr1 target: $CHR1_NAME (label: $LABEL1)"
echo "[INFO] chr2 target: $CHR2_NAME (label: $LABEL2)"
echo "[INFO] PAF: $PAF"
echo "[INFO] Output: $OUT"
echo "[INFO] chr1 threshold: $CHR1_FRAC_THR ; chr2 threshold: $CHR2_FRAC_THR ; dominance: $DOM"

# --------------------
# Step 1: minimap2 (optional)
# --------------------
if [[ "$FILTER_ONLY" -eq 0 ]]; then
  sec_flag="--secondary=no"
  [[ "$SECONDARY" == "yes" ]] && sec_flag="--secondary=yes"
  echo "[INFO] Running minimap2..."
  minimap2 -t "$THREADS" -x "$PRESET" "$sec_flag" "$TARGET" "$QUERY" > "$PAF"
  echo "[INFO] Mapping complete: $PAF"
else
  echo "[INFO] --filter-only enabled; skipping minimap2"
fi

# --------------------
# Step 2: filter/classify from PAF
# --------------------
tmp_hits="$(mktemp)"
tmp_sorted="$(mktemp)"
trap 'rm -f "$tmp_hits" "$tmp_sorted"' EXIT

echo "[INFO] Filtering PAF by target/mapq/alnlen..."
awk -v C1="$CHR1_NAME" -v C2="$CHR2_NAME" -v minMQ="$MIN_MAPQ" -v minAL="$MIN_ALN" '
BEGIN { FS=OFS="\t" }
NF >= 12 {
  t = $6
  if (t != C1 && t != C2) next
  mq = $12 + 0
  al = $11 + 0
  if (mq < minMQ || al < minAL) next

  q = $1
  qlen = $2 + 0
  qs = $3 + 0
  qe = $4 + 0
  if (qe < qs) { tmp = qs; qs = qe; qe = tmp }
  if (qe <= qs) next

  print q, t, qlen, qs, qe, mq, al
}
' "$PAF" > "$tmp_hits"

sort -k1,1 -k2,2 -k4,4n -k5,5n "$tmp_hits" > "$tmp_sorted"

echo "[INFO] Building non-overlapping coverage per scaffold..."
awk -v C1="$CHR1_NAME" -v C2="$CHR2_NAME" -v L1="$LABEL1" -v L2="$LABEL2" -v dom="$DOM" -v c1thr="$CHR1_FRAC_THR" -v c2thr="$CHR2_FRAC_THR" '
BEGIN { FS=OFS="\t" }

function flush_interval(    key, add) {
  if (!have_interval) return
  add = cur_e - cur_s
  if (add > 0) {
    key = cur_q SUBSEP cur_t
    cov[key] += add
  }
}

{
  q = $1
  t = $2
  qlen = $3 + 0
  s = $4 + 0
  e = $5 + 0

  L[q] = qlen
  seen[q] = 1
  hits[q SUBSEP t]++

  key = q SUBSEP t
  if (key != prev_key) {
    flush_interval()
    prev_key = key
    cur_q = q
    cur_t = t
    cur_s = s
    cur_e = e
    have_interval = 1
  } else {
    if (s <= cur_e) {
      if (e > cur_e) cur_e = e
    } else {
      flush_interval()
      cur_s = s
      cur_e = e
      have_interval = 1
    }
  }
}

END {
  flush_interval()

  print "scaffold", "scaffold_len", L1 "_cov_bp", L2 "_cov_bp", L1 "_frac", L2 "_frac", L1 "_hits", L2 "_hits", "call"
  for (q in seen) {
    qlen = L[q] + 0
    c1b = cov[q SUBSEP C1] + 0
    c2b = cov[q SUBSEP C2] + 0
    c1h = hits[q SUBSEP C1] + 0
    c2h = hits[q SUBSEP C2] + 0
    c1f = (qlen > 0 ? c1b / qlen : 0)
    c2f = (qlen > 0 ? c2b / qlen : 0)

    call = "-"
    if (c1b > 0 || c2b > 0) {
      if (c1b > 0 && c2b == 0) {
        if (c1f >= c1thr) call = L1 "_only"
      } else if (c2b > 0 && c1b == 0) {
        if (c2f >= c2thr) call = L2 "_only"
      } else {
        if (c1f >= c1thr && c1f >= c2f * dom) {
          call = L1 "_enriched"
        } else if (c2f >= c2thr && c2f >= c1f * dom) {
          call = L2 "_enriched"
        }
      }
    }

    printf "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d\t%d\t%s\n", q, qlen, c1b, c2b, c1f, c2f, c1h, c2h, call
  }
}
' "$tmp_sorted" > "$OUT"

echo "[INFO] Done. Wrote: $OUT"
