# NuGeneX: Nuclear Gene Extractor and assembler for Fragmented Assemblies

NuGeneX is a feature-rich script for processing, merging, and QC'ing gene contigs against reference sequences. It enables reliable reconstruction of full-length gene sequences from fragmented assemblies and produces detailed quality and summary reports.

## Features

- **Contig Filtering:** Identity and overlap-based filtering and trimming with robust local alignment.
- **Automated Merging:** Sensitive and customizable contig merging with multiple strategies (strict/relaxed).
- **Duplicate Detection:** Evidence-based duplicate gene detection.
- **Comprehensive QC:** Calculates similarity and completeness vs. references, with weighted scoring.
- **Batch Processing:** Processes all genomes in a folder in a single run.
- **Rich Outputs:** Produces FASTA, log, summary, and duplicate evidence files.

---

## Quick Start

You need Python 3 (with [Biopython](https://biopython.org/) and [python-Levenshtein](https://github.com/ztane/python-Levenshtein)):
```sh
pip install biopython python-Levenshtein
```

### **Basic Command**

```sh
python NuGeneX.py -i <input_dir> -r <reference.fasta> -o <output_dir>
```
- `-i` : Folder with input FASTA files (one per genome)
- `-r` : Reference FASTA file (with gene sequences)
- `-o` : Folder for results and logs

### **Complete Example**

```sh
python NuGeneX.py \
  -i input_folder \
  -r resistance_genes.fasta \
  -o results/
```

### **Advanced Options**

- `-s` : Weight for similarity in scoring (default: 6)
- `-c` : Weight for completeness in scoring (default: 4)
- `-m` : Minimum overlap for merging (default: 40)
- `-M` : Minimum overlap for relaxed-merge (default: 5)
- `-l` : Minimum alignment length for duplicate detection (default: 100)
- `-I` : Minimum contig identity (%) (default: 98)

Example:
```sh
python NuGeneX.py \
  -i input_folder \
  -r reference.fasta \
  -o output_folder \
  -s 8 -c 2 -m 30 -M 10 -l 80 -I 95
```

---

## Outputs

- `<genome>_genes_final.fasta`   — Final gene set, best version per gene
- `<genome>_merge_details.log`   — Step-by-step merge logs
- `<genome>_details.txt`         — Duplicate gene evidence
- `<genome>_merged_genes.fasta`  — Merged (multi-contig) gene reconstructions
- `gene_status_summary.csv`      — One-line-per-genome summary table marking Complete/Failed/Duplicated genes

---

**Developed by @zagalom**  
See LICENSE for MIT terms. For questions/issues/collaborations, please open an issue.

## Citation & Contact

Developed by [@zagalom](https://github.com/zagalom).  If you use NuGeneX, please cite the repository.
For questions, bug reports, or collaborations, please open [an issue](https://github.com/zagalom/HapHap/issues).

---

## License

Released under the MIT license. See [LICENSE](LICENSE) for details.
