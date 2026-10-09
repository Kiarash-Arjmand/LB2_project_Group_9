# Laboratory of Bioinformatics II — Group 9

## Project objective

Build and compare protein signal peptide predictors using the von Heijne method and a support vector machine (SVM). Collect annotated eukaryotic protein sequences, reduce sequence redundancy, prepare cross-validation and blind benchmark datasets, analyze the data, implement the methods, and report their performance in a scientific manuscript.

This README describes the current data collection work. Model implementation and evaluation are pending. The workflow follows the supplied course slides (labelled academic year 2025–2026); confirm any updated requirements with the teaching team.

## Current status — 8 October 2026

- [x] Collect UniProt positive and negative candidate annotations.
- [x] Apply annotation filters and retrieve canonical protein sequences.
- [x] Produce the required positive and negative TSV and FASTA files.
- [x] Check IDs, duplicates, lengths, kingdom labels, and negative transmembrane labels.
- [x] Connect to the course VPN.
- [ ] Receive the group's SSH key and server assignment; connect to the course VM.
- [ ] Reduce sequence redundancy with MMseqs2.
- [ ] Select representatives and generate matching TSV files.
- [ ] Make training/benchmark splits and cross-validation folds.
- [ ] Analyze and visualize the datasets.
- [ ] Implement and evaluate von Heijne and SVM predictors.
- [ ] Complete the manuscript and final documentation.

The current datasets are **preliminary, before redundancy reduction**. They must not yet be treated as final training or benchmark sets.

## Local working files

Paths below are relative to the local `LB2_project` folder. The scripts currently sit beside `data/`, at the project root. They have not been verified as uploaded to the shared GitHub repository.

| Path | Purpose |
| --- | --- |
| `filter_positives.py` | Filter positive annotations and prepare positive outputs. |
| `filter_negatives.py` | Filter negative annotations and prepare negative outputs. |
| `data/raw/positive_candidates.tsv` | Original positive candidate annotations. |
| `data/raw/positive_candidates.fasta` | Original positive candidate sequences. |
| `data/raw/positive_candidates_lineage.tsv` | Positive annotations with taxonomic lineage. |
| `data/raw/negative_candidates.tsv` | Original negative candidate annotations. |
| `data/processed/positive.tsv` | Required positive metadata table. |
| `data/processed/positive_filtered.fasta` | Selected canonical positive sequences. |
| `data/processed/negative.tsv` | Required negative metadata table. |
| `data/processed/negative_filtered.fasta` | Selected canonical negative sequences. |
| `data/processed/positive_filtered_annotations.tsv` | Detailed selected positive annotations. |
| `data/processed/positive_excluded_annotations.tsv` | Excluded positive annotations for auditing. |
| `data/processed/negative_filtered_annotations.tsv` | Detailed selected negative annotations. |
| `data/processed/negative_accessions.txt` | Selected negative IDs used for UniProt ID Mapping. |

Preserve the original downloads in `data/raw/`. Store subsequent clustering outputs separately so these preliminary datasets remain available.

## Data source and selection

Source: UniProtKB/Swiss-Prot, retrieved on **6 October 2026**. Counts reflect that retrieval date and may change in later UniProt releases. The UniProt release identifier was not recorded.

Both classes require reviewed eukaryotic proteins (taxonomy ID 2759), no fragments, length at least 40 residues, and protein existence supported at protein level (`existence:1`).

### Positive candidates

UniProt query:

```text
(reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) AND (length:[40 TO *]) AND (existence:1) AND (ft_signal_exp:*)
```

Retain an experimentally supported signal peptide (`ECO:0000269`) starting at residue 1, with a known numeric end position of at least 14. Exclude unknown cleavage positions and features annotated as not cleaved. Store the signal peptide end residue as `cleavage_site`: cleavage occurs after this residue.

### Negative candidates

UniProt query:

```text
(reviewed:true) AND (taxonomy_id:2759) AND (fragment:false) AND (length:[40 TO *]) AND (existence:1) AND NOT (ft_signal:*)
```

Retain proteins with an experimentally supported (`ECO:0000269`) subcellular location in cytosol, nucleus, mitochondrion, plastid, peroxisome, or cell membrane. Evidence must be attached to the relevant location statement. A general Cytoplasm annotation alone does not satisfy the cytosol criterion.

Retrieve selected negative canonical sequences using UniProt ID Mapping: upload `negative_accessions.txt`, map UniProtKB AC/ID to UniProtKB, and download all mapped results as uncompressed canonical FASTA.

### Kingdom assignment

Use taxonomic lineage: Metazoa → `Metazoa`; Fungi → `Fungi`; Viridiplantae → `Plants`; all remaining eukaryotes → `Other`.

## Preliminary dataset counts

| Dataset | Candidates | Retained | Excluded |
| --- | ---: | ---: | ---: |
| Positive | 2,972 | 2,956 | 16 |
| Negative | 66,172 | 16,852 | 49,320 |

| Kingdom | Positive | Negative |
| --- | ---: | ---: |
| Metazoa | 2,440 | 9,985 |
| Fungi | 167 | 3,100 |
| Plants | 312 | 3,513 |
| Other | 37 | 254 |
| Total | 2,956 | 16,852 |

## Required TSV columns

`positive.tsv`:

```text
accession organism kingdom protein_length cleavage_site
```

`negative.tsv`:

```text
accession organism kingdom protein_length tm_helix_first_90
```

Actual files use tabs between columns and include a header row.

The negative flag is `true` when a canonical-sequence TRANSMEM feature has a numeric start within residues **1–90 inclusive** and a note describing a helical segment or signal anchor. Isoform-prefixed coordinates are ignored for canonical FASTA sequences; beta-stranded segments are not counted as helices. An unknown end coordinate does not prevent classification when the start is known. `false` means no qualifying annotated feature was found, rather than experimentally proven absence of a helix.

| Negative transmembrane flag | Proteins |
| --- | ---: |
| `true` | 1,959 |
| `false` | 14,893 |

## Validation completed

- Positive TSV and FASTA each contain 2,956 unique accessions; their ID sets match.
- Negative TSV and FASTA each contain 16,852 unique accessions; their ID sets match.
- No sequence length mismatches were found in either class.
- Positive cleavage positions are at least 14 and less than full protein length.
- The negative table has the required five columns, no blank fields, and kingdom and transmembrane values matching the source annotations and implemented rules.
- No accession is present in both classes. This checks IDs only; similarity between sequences has not yet been assessed.

## Running the existing collection scripts locally

Python 3.13 was used on Windows. The collection scripts use standard-library modules (`csv`, `re`, and `pathlib`). With all expected input files in place, run from the project root:

```bash
python filter_positives.py
python filter_negatives.py
```

These scripts regenerate their output files. The negative FASTA download is a separate UniProt ID Mapping step; `filter_negatives.py` does not download it automatically. The raw UniProt queries and retrieval date above support repeating the collection, although later releases may produce different counts.

## Next preparation stage

According to the supplied preparation slides:

1. Cluster positives and negatives separately with MMseqs2 at 30% sequence identity and 40% coverage, using coverage mode 0 and connected-component clustering (mode 1).
2. Retain one representative per cluster and filter each metadata table to the same representative IDs.
3. Check for sequence similarity across classes before finalizing splits.
4. Assign 80% of each class to training and 20% to the blind benchmark, and generate five stratified cross-validation folds within training. Record the random seed and exact accession assignments.
5. Keep the blind benchmark out of model fitting and hyperparameter selection.

Clustering commands, tool versions, representative counts, split seeds, and output paths will be added after these steps are executed.

## Proposed group responsibilities — not yet agreed

| Person | Proposed responsibility | Deliverables |
| --- | --- | --- |
| Mir Ako Ali | Data collection and annotation QC | Collection scripts, preliminary TSV/FASTA files, documented selection rules. |
| Member 2 — to be assigned | Redundancy reduction, splits, and dataset analysis | Clustering workflow, representative datasets, shared splits/folds, descriptive figures. |
| Member 3 — to be assigned | von Heijne predictor | Training/scoring code, validation-based threshold selection, evaluation results. |
| Member 4 — to be assigned | SVM predictor | Feature extraction, hyperparameter selection, training and evaluation results. |
| All four members | Integration, interpretation, manuscript, and oral preparation | Common evaluation protocol, reviewed code, manuscript, understanding of the full workflow. |

Agree these roles together before treating them as assignments. Both predictors must use the same final datasets, folds, and benchmark accession lists to support a fair comparison.

## Access and confidentiality

The course VM is reached through the course VPN and authenticated using the professor-provided group SSH private key. VPN connection is working; SSH credentials and server assignment are still pending.

Do not commit VPN profiles (`*.ovpn`), SSH private keys, passwords, or other credentials to GitHub. Share the project's data and code through the group repository once its current structure and existing files have been checked. This README does not indicate that a GitHub upload has already occurred.
