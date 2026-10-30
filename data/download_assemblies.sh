#!/usr/bin/env bash
# Download the 10 Salmonella enterica genome assemblies used in this project.
# Requires the NCBI Datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download/
# Usage: bash data/download_assemblies.sh   (run from the repo root)
# Output: data/assemblies/<LABEL>.fna  (one FASTA per genome)
set -euo pipefail

OUTDIR="data/assemblies"
mkdir -p "$OUTDIR"

# accession, short label, serovar, strain
GENOMES=(
  "GCF_000006945.2 Typhimurium_LT2"
  "GCF_000210855.2 Typhimurium_SL1344"
  "GCF_008244845.1 Typhimurium_14028S"
  "GCF_000009505.1 Enteritidis_P125109"
  "GCF_000195995.1 Typhi_CT18"
  "GCF_000007545.1 Typhi_Ty2"
  "GCF_000020705.1 Heidelberg_SL476"
  "GCF_000016045.1 Newport_SL254"
  "GCF_000008105.1 Choleraesuis_SC-B67"
  "GCA_010634545.1 ParatyphiA_ATCC9150"
)

for entry in "${GENOMES[@]}"; do
  acc="${entry%% *}"
  label="${entry##* }"
  echo "Downloading $acc ($label)..."
  datasets download genome accession "$acc" \
    --include genome \
    --filename "/tmp/${acc}.zip"
  unzip -o -q "/tmp/${acc}.zip" -d "/tmp/${acc}_unzipped"
  fna=$(find "/tmp/${acc}_unzipped" -name "*.fna" | head -n 1)
  cp "$fna" "$OUTDIR/${label}.fna"
  rm -rf "/tmp/${acc}.zip" "/tmp/${acc}_unzipped"
  echo "  saved to $OUTDIR/${label}.fna"
done

echo "Done. $(ls "$OUTDIR"/*.fna | wc -l) assemblies in $OUTDIR/"
