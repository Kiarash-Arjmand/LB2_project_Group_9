# Laboratory of Bioinformatics II — Group 9

We are comparing the von Heijne method and a support vector machine for predicting secretory signal peptides in eukaryotic proteins.

## Current progress

The datasets contain **2,956 positives** and **20,974 negatives**, before clustering. The negative collection was corrected on 12 October 2026. Sequence clustering, dataset analysis, training and evaluation are still to be done.

## Data collection

Both classes contain reviewed eukaryotic proteins with evidence at protein level, no fragment annotation and a length of at least 40 residues.

Positives were collected on 6 October 2026. The positive filter keeps experimentally supported signal peptides starting at residue 1, with a known cleavage position and a length of at least 14 residues. Non-cleaved features are excluded. The cleavage position is the last signal peptide residue.

Negatives were recollected on 12 October 2026 from **UniProt release 2026_03**. The search excludes signal peptides at any evidence level and uses UniProt's experimental subcellular-location search for cytosol, nucleus, mitochondrion, plastid, peroxisome or cell membrane. The exact query, release and file checksums are recorded in `data/raw/negative_collection.json`.

The earlier negative search omitted the location restriction and returned 66,172 candidates. Its local annotation filter retained 16,852. The corrected search returns 20,974; the replacement script checks and formats those results without applying the old location regex again. All 16,852 previous negatives remain in the corrected set, and 4,122 have been added.

## Files

Paths are relative to the project root. Scripts sit beside `data/`.

| Path | Contents |
| --- | --- |
| `filter_positives.py` | Positive filtering script. |
| `collect_negative_candidates.py` | Download the course-specific negative collection and prepare its outputs. |
| `filter_negatives.py` | Check the recorded negative inputs and regenerate the processed files. |
| `data/raw/positive_candidates.tsv` | Positive candidate annotations. |
| `data/raw/positive_candidates.fasta` | Positive candidate sequences. |
| `data/raw/positive_candidates_lineage.tsv` | Positive annotations with taxonomic lineage. |
| `data/raw/negative_candidates.tsv` | Corrected negative annotations, including protein existence and fragment status. |
| `data/raw/negative_candidates.fasta` | Corresponding canonical negative sequences. |
| `data/raw/negative_collection.json` | Negative collection query and provenance. |
| `data/processed/positive.tsv` | Positive metadata. |
| `data/processed/positive_filtered.fasta` | Retained positive sequences. |
| `data/processed/positive_filtered_annotations.tsv` | Retained positive annotations. |
| `data/processed/positive_excluded_annotations.tsv` | Excluded positive annotations. |
| `data/processed/negative.tsv` | Negative metadata. |
| `data/processed/negative_filtered.fasta` | Negative sequences. |
| `data/processed/negative_filtered_annotations.tsv` | Detailed negative annotations. |
| `data/processed/negative_accessions.txt` | Negative accession list. |
| `data/audit/` | Counts, checks and added accessions from the negative correction. |

## Dataset counts

| Dataset | Candidates | Retained | Excluded |
| --- | ---: | ---: | ---: |
| Positive | 2,972 | 2,956 | 16 |
| Negative | 20,974 | 20,974 | 0 |

| Kingdom | Positive | Negative |
| --- | ---: | ---: |
| Metazoa | 2,440 | 12,634 |
| Fungi | 167 | 3,814 |
| Plants | 312 | 4,132 |
| Other | 37 | 394 |
| Total | 2,956 | 20,974 |

## Metadata

The metadata tables are TSV files with headers.

| File | Columns |
| --- | --- |
| `positive.tsv` | `accession`, `organism`, `kingdom`, `protein_length`, `cleavage_site` |
| `negative.tsv` | `accession`, `organism`, `kingdom`, `protein_length`, `tm_helix_first_90` |

Kingdoms are assigned from taxonomic lineage: Metazoa, Fungi, Plants (Viridiplantae), or Other.

For negatives, `tm_helix_first_90` is `true` when a canonical helical transmembrane segment or signal anchor starts within residues 1–90. Isoform-specific coordinates and beta-stranded segments are excluded. There are 2,533 `true` values and 18,441 `false` values.

The corrected negative TSV and FASTA have matching unique accessions and sequence lengths. There is no accession overlap with the current positive table. Existing negatives retain the same metadata and sequences.

## Reproducing the negative files

Python's standard library is sufficient. From the project root, run:

```bash
python filter_negatives.py
```

This uses the included raw files and collection record, checks them and regenerates the four processed negative files. Existing processed files are backed up under `data/archive/`.

To repeat the UniProt download instead, run:

```bash
python collect_negative_candidates.py
```

This requires internet access. It downloads all search results, retrieves the corresponding canonical sequences and regenerates the processed files. Existing raw and processed negatives are backed up. Later UniProt releases may return different counts.

## Next step

Cluster the positive and negative sequences separately with MMseqs2, select representatives, then create the training and benchmark sets and cross-validation folds. The current counts are collection counts; they will decrease after redundancy reduction.
