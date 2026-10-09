# Laboratory of Bioinformatics II — Group 9

This project compares two methods for predicting secretory signal peptides in eukaryotic proteins: the von Heijne method and a support vector machine (SVM).

## Progress

Data collection and annotation filtering are complete. The filtering scripts and processed datasets have been uploaded to the group repository.

The current dataset contains **2,956 positive proteins** and **16,852 negative proteins**. These counts are before sequence clustering. Redundancy reduction, training/test splits, data analysis, and predictor development are still pending.

## Files

All paths are relative to the repository root. The two filtering scripts are at the root, beside `data/`.

| File or folder | Contents |
| --- | --- |
| `filter_positives.py` | Positive annotation filtering and output preparation. |
| `filter_negatives.py` | Negative annotation filtering and output preparation. |
| `data/raw/positive_candidates.tsv` | Original positive candidate annotations. |
| `data/raw/positive_candidates.fasta` | Original positive candidate sequences. |
| `data/raw/positive_candidates_lineage.tsv` | Positive annotations with taxonomic lineage. |
| `data/raw/negative_candidates.tsv` | Original negative candidate annotations. |
| `data/processed/positive.tsv` | Metadata for retained positive proteins. |
| `data/processed/positive_filtered.fasta` | Full canonical sequences of retained positive proteins. |
| `data/processed/negative.tsv` | Metadata for retained negative proteins. |
| `data/processed/negative_filtered.fasta` | Full canonical sequences of retained negative proteins. |
| `data/processed/positive_filtered_annotations.tsv` | Detailed annotations for retained positives. |
| `data/processed/positive_excluded_annotations.tsv` | Positive annotations excluded during filtering. |
| `data/processed/negative_filtered_annotations.tsv` | Detailed annotations for retained negatives. |
| `data/processed/negative_accessions.txt` | Negative accession list used to retrieve sequences. |
| `data/splits/` | Reserved for training, benchmark, and fold assignments; not generated yet. |
| `src/`, `notebooks/`, `figures/` | Reserved for subsequent code, analyses, and figures. |

## Data collection

The candidate data were downloaded from **UniProtKB/Swiss-Prot on 6 October 2026**. The release identifier was not recorded, so repeating the searches on a later release may give different counts.

Both searches select reviewed eukaryotic proteins (taxonomy ID 2759), exclude fragments, require a sequence length of at least 40 residues, and require evidence at protein level (`existence:1`).

### Positive proteins

Search query:

```text
(reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) AND (length:[40 TO *]) AND (existence:1) AND (ft_signal_exp:*)
```

The positive filter retains signal peptides with experimental evidence (`ECO:0000269`) that start at residue 1 and have a known end position of at least 14. Features with unknown cleavage positions or a “not cleaved” annotation are excluded.

The `cleavage_site` value is the last residue of the signal peptide: cleavage occurs immediately after it.

### Negative proteins

Search query:

```text
(reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) AND (length:[40 TO *]) AND (existence:1) AND NOT (ft_signal:*)
```

The negative filter retains proteins with experimentally supported localization (`ECO:0000269`) in the cytosol, nucleus, mitochondrion, plastid, peroxisome, or cell membrane. The evidence must apply to the relevant location statement. “Cytoplasm” alone does not meet the cytosol criterion.

Canonical negative sequences were retrieved through UniProt ID Mapping using `negative_accessions.txt`.

### Taxonomic groups

Taxonomic lineage is mapped to four groups: Metazoa, Fungi, Plants (Viridiplantae), and Other.

## Filtering results

| Dataset | Candidates | Retained | Excluded |
| --- | ---: | ---: | ---: |
| Positive | 2,972 | 2,956 | 16 |
| Negative | 66,172 | 16,852 | 49,320 |

| Group | Positive | Negative |
| --- | ---: | ---: |
| Metazoa | 2,440 | 9,985 |
| Fungi | 167 | 3,100 |
| Plants | 312 | 3,513 |
| Other | 37 | 254 |
| Total | 2,956 | 16,852 |

## Output format

Both metadata files are tab-separated and include a header row.

| File | Columns, in order |
| --- | --- |
| `positive.tsv` | `accession`, `organism`, `kingdom`, `protein_length`, `cleavage_site` |
| `negative.tsv` | `accession`, `organism`, `kingdom`, `protein_length`, `tm_helix_first_90` |

In the negative table, `tm_helix_first_90` is `true` when an annotated helical transmembrane segment or signal anchor starts within residues 1–90, inclusive. The filter uses canonical-sequence coordinates, ignores isoform-prefixed features, and excludes beta-stranded segments. A known start can be used even when the end position is unknown.

There are **1,959** negatives flagged `true` and **14,893** flagged `false`. A `false` value means that no qualifying feature was found in the annotations.

## Checks completed

- Each metadata table and its corresponding FASTA contain the same unique accession IDs.
- Sequence lengths match the metadata in both classes.
- Positive cleavage positions are at least 14 and below the full sequence length.
- The negative table has five columns with no blank fields; taxonomic groups and transmembrane flags match the filtering rules.
- No accession appears in both classes.

Sequence similarity has not yet been assessed. These files are the inputs to redundancy reduction, not the final training and benchmark datasets.

## Running the filters

The filtering scripts use Python's standard library; no additional packages are needed for this step. Place the input files in `data/raw/`, then run these commands from the repository root:

```bash
python filter_positives.py
python filter_negatives.py
```

Running the scripts regenerates their output files. Negative sequence retrieval is a separate step: map the IDs in `data/processed/negative_accessions.txt` from UniProtKB AC/ID to UniProtKB, download canonical FASTA sequences, and save them as `data/processed/negative_filtered.fasta`.

## Next steps

Reduce sequence redundancy with MMseqs2, select representative proteins, and create the training/benchmark split and cross-validation folds. Dataset analysis and implementation of the von Heijne and SVM methods will follow. Commands, parameters, and results will be documented as each stage is completed.
