# NuGeneX
Nuclear Gene Extractor for Fragmented Assemblies. NuGeneX is a Python-based tool designed to extract and assemble full-length nuclear gene sequences from fragmented assemblies. It maps each target locus to the reference supplied and reconstructs contiguous sequences for downstream analyses.
Summary:

    Purpose: The script automates the process of refining gene assemblies (from FASTA files)—detecting duplicates, trimming, merging contigs, assessing sequence quality, and generating summary reports, with numerous customizable parameters for thresholds and scoring.
    How it Works:
        Input: Takes a directory of input FASTA files (assembled gene contigs), a reference gene FASTA, and outputs results to a specified directory.
        Processing Includes:
            Trimming contigs to match reference gene boundaries using local alignment.
            Splitting contigs with ambiguous (N) regions.
            Filtering contigs by minimum identity threshold.
            Detecting potential gene duplications based on sequence similarity and overlap.
            Attempting to merge overlapping contigs using several strategies (strict, standard, relaxed).
            Scoring each sequence (or merged chain) based on similarity and completeness regarding the reference.
        Outputs:
            Best reconstructed gene version(s) for each input genome—either a single contig or a merged sequence—with annotated quality.
            Logs detailing merge and duplication detection steps.
            A summary CSV file reporting status ("Complete", "Failed", "Duplicated", etc.) for each gene across all processed genomes.

Key Features:

    Automated handling of ambiguous regions, duplicate detection, and multi-round progressive merging of contigs.
    Highly configurable through command-line arguments: scoring weights, overlap size thresholds, identity thresholds, etc.
    Generates both detailed logs for troubleshooting and high-level CSV summaries for comparing gene assembly completeness/quality across samples.

Typical Use Case:
Researchers who assemble draft genomes and want to reliably reconstruct AMR gene sequences (or similar genes) by merging fragmented contigs, detecting possible duplications, and scoring overall quality relative to reference genes.
