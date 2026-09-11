# Merge Metschnikowia OTU features (bacterial/fungal feature table, taxonomy
# and representative sequences) according to the phylogeny-informed merge
# map produced by metschnikowia_tree_merge.R.

import os
import pandas as pd
from Bio import SeqIO

# --- Import data ---------------------------------------------------------
# metschnikowia_merge_map.tsv: feature_id -> merged_id, from metschnikowia_tree_merge.R
# otutable.txt: QIIME2-exported feature table (TSV, one comment line to skip)
# taxonomy.tsv: QIIME2 taxonomy assignments
# dna-sequences.fasta: representative sequences, one per feature
data_dir = "data"
output_dir = "output"

merge_map_path = os.path.join(data_dir, "metschnikowia_merge_map.tsv")
table_path = os.path.join(data_dir, "otutable.txt")
taxo_path = os.path.join(data_dir, "taxonomy.tsv")
seqs_path = os.path.join(data_dir, "dna-sequences.fasta")

# --- 1. Load the merge map ------------------------------------------------
merge_map = pd.read_csv(merge_map_path, sep="\t")
merge_dict = dict(zip(merge_map["feature_id"], merge_map["merged_id"]))

print(f"Merge map: {len(merge_map)} entries")
print(f"  OTUs before merge: {merge_map['feature_id'].nunique()}")
print(f"  OTUs after merge: {merge_map['merged_id'].nunique()}")

# --- 2. Merge the feature table --------------------------------------------
df = pd.read_csv(table_path, sep="\t", skiprows=1, index_col=0)
print(f"\nFeature table loaded: {df.shape[0]} features x {df.shape[1]} samples")

# rename features per the merge map; non-Metschnikowia features are unchanged
df.index = df.index.map(lambda x: merge_dict.get(x, x))
df_merged = df.groupby(df.index).sum()

print(f"Feature table after merge: {df_merged.shape[0]} features x {df_merged.shape[1]} samples")
print(f"Features reduced: {df.shape[0]} -> {df_merged.shape[0]} (-{df.shape[0] - df_merged.shape[0]})")

output_table = os.path.join(output_dir, "otutable_merged.txt")
df_merged.to_csv(output_table, sep="\t")
print(f"Exported: {output_table}")

# --- 3. Filter the taxonomy to the merged feature set -----------------------
taxo = pd.read_csv(taxo_path, sep="\t")
print(f"\nTaxonomy loaded: {len(taxo)} entries")

merged_ids = set(df_merged.index)
taxo_filtered = taxo[taxo.iloc[:, 0].isin(merged_ids)]
print(f"Taxonomy after filtering: {len(taxo_filtered)} entries")

output_taxo = os.path.join(output_dir, "taxonomy_merged.tsv")
taxo_filtered.to_csv(output_taxo, sep="\t", index=False)
print(f"Exported: {output_taxo}")

# --- 4. Filter representative sequences to the merged feature set -----------
seqs = SeqIO.to_dict(SeqIO.parse(seqs_path, "fasta"))
print(f"\nSequences loaded: {len(seqs)}")

kept_seqs = [seqs[sid] for sid in merged_ids if sid in seqs]
print(f"Sequences after filtering: {len(kept_seqs)}")

output_seqs = os.path.join(output_dir, "rep-seqs-merged.fasta")
SeqIO.write(kept_seqs, output_seqs, "fasta")
print(f"Exported: {output_seqs}")

# --- 5. Summary --------------------------------------------------------------
print("\n=== SUMMARY ===")
print(f"Feature table: {df.shape[0]} -> {df_merged.shape[0]} features")
print(f"Taxonomy: {len(taxo)} -> {len(taxo_filtered)} entries")
print(f"Sequences: {len(seqs)} -> {len(kept_seqs)}")
print(f"\nFiles written to {output_dir}:")
print("  - otutable_merged.txt")
print("  - taxonomy_merged.tsv")
print("  - rep-seqs-merged.fasta")
print("\nTo re-import into QIIME2:")
print("  qiime tools import --input-path otutable_merged.txt --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path table-merged.qza")
print("  qiime tools import --input-path taxonomy_merged.tsv --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --output-path taxonomy-merged.qza")
print("  qiime tools import --input-path rep-seqs-merged.fasta --type 'FeatureData[Sequence]' --output-path rep-seqs-merged.qza")