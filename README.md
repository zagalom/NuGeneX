#  NuGeneX: Nuclear Gene Extractor and assembler for Fragmented Assemblies

**NuGeneX** is a four-step automated pipeline for the reconstruction and evaluation of nuclear genes from draft assemblies. In its development used *Candida albicans* (and other yeasts) as a models.  
It integrates BLAST-based gene detection, contig extraction, overlap assembly, and sequence evaluation into a fully automated workflow optimized for large datasets.

---

##  Overview

NuGeneX automates the retrieval of gene sequences from multiple assemblies using a reference CDS file.  
It performs similarity-based searches with BLAST, extracts the matching regions, merges overlapping fragments, and evaluates each reconstructed gene based on similarity and completeness.

### Main Features
- 🔹 Automated BLAST database creation and querying  
- 🔹 Smart contig extraction with identity and length thresholds  
- 🔹 Consolidation of fragmented hits on the same scaffold  
- 🔹 Local alignment-based trimming for reference-bound sequences  
- 🔹 Recursive overlap merging and scoring  
- 🔹 Parallel-ready design for large genomic datasets  

---

##  Pipeline Structure

| Step | Description | Main Function |
|------|--------------|----------------|
| **1** | Create nucleotide BLAST databases for each assembly | `step1_create_blast_databases()` |
| **2** | Run BLAST searches against the reference gene set | `step2_run_blast()` |
| **3** | Extract candidate contigs from scaffolds based on BLAST hits | `step3_extract_contigs()` |
| **4** | Merge overlapping contigs, trim, score, and select best assemblies | `step4_assemble_and_evaluate()` |

---

##  Example Usage

```bash
# STEP 1: Create BLAST databases
python NuGeneX.py step1 -i assemblies/ -o blast_dbs/

# STEP 2: Run BLAST using the reference CDS file as query
python NuGeneX.py step2 -d blast_dbs/ -q reference_CDS.fasta -o blast_results/

# STEP 3: Extract candidate contigs from scaffolds
python NuGeneX.py step3 -r reference_CDS.fasta -b blast_results/ -s scaffolds/ -o extracted_genes/

# STEP 4: Assemble and evaluate genes
python NuGeneX.py step4 -i extracted_genes/ -r reference_CDS.fasta -o final_assemblies/
```

Each step can be run independently — perfect for distributed or checkpointed execution in large projects.

---

##  Input Requirements

| File Type | Description |
|------------|--------------|
| **Reference CDS file** | FASTA file containing all reference gene coding sequences |
| **Scaffold FASTAs** | Genome assembly files (e.g. `sample1.fasta`, `sample2.fasta`, ...) |
| **BLAST results** | Automatically generated in Step 2, used by Step 3 |

---

##  Output Files

| Output Directory | Description |
|------------------|--------------|
| **blast_dbs/** | Local nucleotide BLAST databases |
| **blast_results/** | Raw BLAST results per genome |
| **extracted_genes/** | Per-genome FASTA files containing extracted contigs |
| **final_assemblies/** | Best reconstructed sequences per gene and genome |

Each record header encodes key metadata for downstream analysis:  
```
GeneID|ScaffoldID|RefStart-RefEnd|HitStart-HitEnd|Identity|Length|[Consolidated]
```

---

##  Scoring System

NuGeneX evaluates each assembled gene using a weighted score:

\[
\text{Score} = \frac{(Similarity × 6) + (Completeness × 4)}{10}
\]

- **Similarity** → based on modified Levenshtein distance  
- **Completeness** → relative to the full reference gene length  
- **Penalty** for over-completeness (>100%) is applied automatically  

---

##  Dependencies

- Python ≥ 3.9  
- [Biopython](https://biopython.org)  
- [Levenshtein](https://pypi.org/project/python-Levenshtein/)  
- NCBI BLAST+ tools (`makeblastdb`, `blastn`)  

Install dependencies:

```bash
pip install biopython python-Levenshtein
```

---

## Citation & Context

This tool was designed to support large-scale reconstruction of nuclear gene sequences of multiple isolates, particularly for studies on population genomics.
If you use NuGeneX, please cite the repository.
> 
> (Manuscript in preparation).

---

##  Author

Developed by [@zagalom](https://github.com/zagalom).  If you use NuGeneX, please cite the repository.
For questions, bug reports, or collaborations, please open [an issue](https://github.com/zagalom/HapHap/issues).

---

##  License


Released under the MIT license. See [LICENSE](LICENSE) for details.

---

