import os
import subprocess
import argparse
import re
import sys
import shutil
import platform
from functools import lru_cache

# Check for required Python packages and provide helpful instructions if missing
_missing_packages = []
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import PairwiseAligner
except ModuleNotFoundError:
    _missing_packages.append("biopython (pip install biopython)")

try:
    from Levenshtein import distance as levenshtein_distance
except ModuleNotFoundError:
    _missing_packages.append("python-Levenshtein (pip install python-Levenshtein)")

if _missing_packages:
    print("❌ Missing required Python packages: " + ", ".join(_missing_packages))
    print("\nPlease install them. Example (Windows):")
    print("  python -m pip install --upgrade pip")
    print("  python -m pip install biopython python-Levenshtein")
    print("\nOr using a virtual environment:")
    print("  python -m venv venv")
    if platform.system() == "Windows":
        print("  .\\venv\\Scripts\\activate")
    else:
        print("  source venv/bin/activate")
    print("  pip install biopython python-Levenshtein")

    # Check for BLAST+ tools and provide a hint if missing
    missing_tools = []
    for tool in ("makeblastdb", "blastn"):
        if shutil.which(tool) is None:
            missing_tools.append(tool)

    if missing_tools:
        print("\nAlso, the following external BLAST+ tools were not found on your PATH: " + ", ".join(missing_tools))
        print("If you need to run BLAST from this script, download BLAST+ from:")
        print("  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
        print("Then add the directory containing makeblastdb and blastn to your PATH.")
    sys.exit(1)

# ==============================================================================
# --- STEP 1: CREATE BLAST DATABASES ---
# ==============================================================================

def step1_create_blast_databases(input_dir: str, output_dir: str):
    """
    Step 1: Creates a BLAST nucleotide database for all FASTA files
            in the input directory.
    """
    print("🚀 Starting Step 1: Creating BLAST Databases...")
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]

    for fasta_file in fasta_files:
        fasta_path = os.path.join(input_dir, fasta_file)
        db_name = os.path.join(output_dir, fasta_file.replace(".fasta", ""))
        print(f"  -> Creating database for: {fasta_file}...")

        try:
            subprocess.run([
                "makeblastdb",
                "-in", fasta_path,
                "-dbtype", "nucl",
                "-out", db_name
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Error creating DB for {fasta_file}: {e}")
            continue
        print(f"  ✅ Database created: {os.path.basename(db_name)}")

    print("✅ Step 1 completed: All necessary BLAST databases created.")

# ==============================================================================
# --- STEP 2: RUN BLAST SEARCH ---
# ==============================================================================

def step2_run_blast(db_dir: str, query_file: str, output_dir: str):
    """
    Step 2: Runs blastn using the query file against all databases in db_dir.
    """
    print("\n🚀 Starting Step 2: Running BLAST Searches...")
    os.makedirs(output_dir, exist_ok=True)
    db_files = [f for f in os.listdir(db_dir) if f.endswith(".nhr")]

    for db_file in db_files:
        db_name = os.path.join(db_dir, db_file.replace(".nhr", ""))
        # This output name must match what Step 3 expects
        output_file = os.path.join(output_dir, f"{os.path.basename(db_name)}_results.txt")

        print(f"  -> Running BLAST on database: {os.path.basename(db_name)}...")

        try:
            subprocess.run([
                "blastn",
                "-query", query_file,
                "-db", db_name,
                "-out", output_file,
                "-outfmt", "6",
                "-evalue", "1e-5",
                "-max_target_seqs", "10"
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Error running BLAST on {os.path.basename(db_name)}: {e}")
            continue
        print(f"  ✅ Results saved to: {os.path.basename(output_file)}")

    print("✅ Step 2 completed: All BLAST searches executed.")

# ==============================================================================
# --- STEP 3: EXTRACT CONTIG SEQUENCES ---
# ==============================================================================

def read_reference_lengths(reference_file):
    """Reads gene lengths from the reference FASTA."""
    reference_lengths = {}
    try:
        for record in SeqIO.parse(reference_file, "fasta"):
            reference_lengths[record.id] = len(record.seq)
    except FileNotFoundError:
        print(f"Error: Reference file not found at {reference_file}")
        raise
    return reference_lengths

def group_hits_by_node(filtered_hits):
    """Groups and consolidates hits on the same scaffold (NODE) for a single query."""
    node_groups = {}

    for hit in filtered_hits:
        key = (hit['query_id'], hit['scaffold_id'])
        if key not in node_groups:
            node_groups[key] = []
        node_groups[key].append(hit)

    consolidated_hits = []
    for (query_id, scaffold_id), hits in node_groups.items():
        if len(hits) >= 2:
            strand = '-' if hits[0]['hit_start'] > hits[0]['hit_end'] else '+'
            best_pident = max(hit['pident'] for hit in hits)

            all_ref_starts = [hit['ref_start'] for hit in hits]
            all_ref_ends = [hit['ref_end'] for hit in hits]
            all_hit_starts = [hit['hit_start'] for hit in hits]
            all_hit_ends = [hit['hit_end'] for hit in hits]

            min_ref = min(min(all_ref_starts), min(all_ref_ends))
            max_ref = max(max(all_ref_starts), max(all_ref_ends))
            min_hit = min(min(all_hit_starts), min(all_hit_ends))
            max_hit = max(max(all_hit_starts), max(all_hit_ends))

            consolidated_hits.append({
                'query_id': query_id,
                'scaffold_id': scaffold_id,
                'pident': best_pident,
                'align_length': max_hit - min_hit,
                'ref_start': min_ref,
                'ref_end': max_ref,
                'hit_start': min_hit if strand == '+' else max_hit,
                'hit_end': max_hit if strand == '+' else min_hit,
                'is_consolidated': True
            })
        else:
            consolidated_hits.extend(hits)

    return consolidated_hits

def extract_sequences_from_blast(blast_file, fasta_file, output_file, reference_lengths, min_pident=70.0, min_length=50):
    """Extracts sequences from scaffolds based on BLAST hits and writes them to FASTA."""
    try:
        scaffold_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    except FileNotFoundError:
        print(f"  Warning: Scaffold FASTA not found at {fasta_file}, skipping extraction.")
        return

    blast_hits = []
    try:
        with open(blast_file, "r") as bf:
            for line in bf:
                cols = line.strip().split("\t")
                if len(cols) < 10: continue
                query_id, scaffold_id = cols[0], cols[1]
                pident = float(cols[2])
                align_length = int(cols[3])
                ref_start, ref_end = int(cols[6]), int(cols[7])
                hit_start, hit_end = int(cols[8]), int(cols[9])

                if pident < min_pident or align_length < min_length:
                    continue

                blast_hits.append({
                    'query_id': query_id,
                    'scaffold_id': scaffold_id,
                    'pident': pident,
                    'align_length': align_length,
                    'ref_start': min(ref_start, ref_end),
                    'ref_end': max(ref_start, ref_end),
                    'hit_start': hit_start,
                    'hit_end': hit_end
                })
    except FileNotFoundError:
        print(f"  Warning: BLAST results file not found at {blast_file}, skipping extraction.")
        return

    blast_hits.sort(key=lambda x: (-x['pident'], -x['align_length']))

    filtered_hits = []
    seen_regions = {}
    for hit in blast_hits:
        # Use query_id, scaffold_id, and ref_start as a key to group overlapping fragments
        gene_ref_key = (hit['query_id'], hit['scaffold_id'], hit['ref_start'])

        # This logic is simplified: it keeps hits if they are not in a seen region
        # or have a higher pident than a previously seen hit in that *exact* region.
        # A more complex logic would merge overlapping regions.
        if gene_ref_key not in seen_regions or hit['pident'] > seen_regions[gene_ref_key]['pident']:
             if gene_ref_key in seen_regions:
                 filtered_hits = [h for h in filtered_hits if (h['query_id'], h['scaffold_id'], h['ref_start']) != gene_ref_key]

             filtered_hits.append(hit)
             seen_regions[gene_ref_key] = hit

    filtered_hits.sort(key=lambda x: (x['query_id'], x['ref_start']))
    final_hits = group_hits_by_node(filtered_hits)

    extracted_sequences = []
    for hit in final_hits:
        query_id = hit['query_id']
        scaffold_id = hit['scaffold_id']
        pident = hit['pident']
        ref_start = hit['ref_start']
        ref_end = hit['ref_end']
        hit_start = hit['hit_start']
        hit_end = hit['hit_end']

        gene_length = reference_lengths.get(query_id)
        if not gene_length: continue

        if scaffold_id in scaffold_sequences:
            scaffold_length = len(scaffold_sequences[scaffold_id])
            strand = '-' if hit_start > hit_end else '+'
            start, end = min(hit_start, hit_end), max(hit_start, hit_end)

            if strand == '+':
                new_start = max(1, start - (ref_start - 1))
                new_end = min(scaffold_length, end + (gene_length - ref_end))
            else:
                new_start = max(1, start - (gene_length - ref_end))
                new_end = min(scaffold_length, end + (ref_start - 1))

            seq = scaffold_sequences[scaffold_id].seq[new_start-1:new_end]
            if strand == '-':
                seq = seq.reverse_complement()

            current_len = len(seq)
            # This is the new header format that Step 4 must parse
            header_parts = [
                query_id,
                scaffold_id,
                f"{ref_start}-{ref_end}",
                f"{new_start}-{new_end}",
                f"{pident:.2f}",
                str(current_len)
            ]

            if hit.get('is_consolidated'):
                header_parts.append('Consolidated')

            header = "|".join(header_parts)
            if header and seq:
                extracted_sequences.append(SeqRecord(seq, id=header, description=""))

    extracted_sequences.sort(key=lambda x: (x.id.split("|")[0], float(x.id.split("|")[4]), len(x.seq)), reverse=True)

    try:
        SeqIO.write(extracted_sequences, output_file, "fasta")
    except Exception as e:
        print(f"  ❌ Error writing output file {output_file}: {e}")

def step3_extract_contigs(blast_results_dir: str, scaffold_fasta_dir: str, output_fasta_dir: str, reference_cds_file: str):
    """
    Step 3: Extracts candidate gene sequences from scaffolds based on BLAST results.
    """
    print("\n🚀 Starting Step 3: Extracting Contigs from Scaffolds...")
    os.makedirs(output_fasta_dir, exist_ok=True)

    try:
        reference_lengths = read_reference_lengths(reference_cds_file)
    except FileNotFoundError:
        print("❌ Cannot proceed with Step 3: Reference CDS file not found.")
        return

    # Step 3 must look for the file names Step 2 created.
    for blast_file in os.listdir(blast_results_dir):
        if blast_file.endswith("_results.txt"):
            base_name = blast_file.replace("_results.txt", "")
            scaffold_fasta = os.path.join(scaffold_fasta_dir, f"{base_name}.fasta")
            # This output name must match what Step 4 expects
            output_fasta = os.path.join(output_fasta_dir, f"{base_name}_genes.fasta")

            print(f"  -> Processing BLAST results for: {base_name}")
            extract_sequences_from_blast(
                blast_file=os.path.join(blast_results_dir, blast_file),
                fasta_file=scaffold_fasta,
                output_file=output_fasta,
                reference_lengths=reference_lengths
            )
            print(f"  ✅ Contig file created: {os.path.basename(output_fasta)}")

    print("✅ Step 3 completed: All contigs extracted and saved.")

# ==============================================================================
# --- STEP 4: ASSEMBLE & EVALUATE GENES ---
# (This section now contains ALL functions from your original script)
# ==============================================================================

@lru_cache(maxsize=None)
def calculate_similarity(reference_seq, merged_seq, ref_start, ref_end):
    """Calculate similarity using modified Levenshtein distance"""
    if not reference_seq or not merged_seq:
        return 0.0

    # Get the relevant part of reference sequence
    reference_region = str(reference_seq[ref_start:ref_end])
    merged_seq = str(merged_seq)

    if not reference_region or not merged_seq:
        return 0.0

    # Calculate Levenshtein similarity
    max_length = max(len(reference_region), len(merged_seq))
    if max_length == 0:
        return 100.0 # Both are empty, perfect match

    distance = levenshtein_distance(reference_region, merged_seq)
    similarity = (1 - (distance / max_length)) * 100

    return similarity

def calculate_score(similarity, completeness, score_weights=(6, 4, 10)):
    """Calculate merge score with penalty for over-completeness (more than 100%)"""
    if completeness > 100:
        penalty = (completeness - 100) * 1 # Penalty of 1 per percentage point
        adjusted_completeness = 100 - penalty
    else:
        adjusted_completeness = completeness

    # Weight similarity more heavily than completeness
    sim_weight, comp_weight, total = score_weights
    return (similarity * sim_weight + adjusted_completeness * comp_weight) / total

def split_sequence_on_n(record, reference_seq):
    """Split sequences containing N regions and trim resulting fragments"""
    seq_str = str(record.seq)
    n_regions = list(re.finditer('N{5,}', seq_str))

    if not n_regions:
        # No Ns found, just trim the original record
        trimmed = trim_to_reference_bounds(record, reference_seq)
        return [trimmed] if trimmed else []

    sequences = []
    parts = record.id.split("|")
    gene_name = parts[0]
    node_base = parts[1]

    # Get original positions and identity
    orig_ref_start = int(parts[2].split("-")[0])
    orig_ref_end = int(parts[2].split("-")[1])
    orig_query_start = int(parts[3].split("-")[0])
    orig_query_end = int(parts[3].split("-")[1])
    identity = parts[4]

    # Process first fragment
    first_fragment_str = seq_str[:n_regions[0].start()]
    if first_fragment_str:
        first_fragment_length = len(first_fragment_str)
        new_ref_end = orig_ref_start + first_fragment_length - 1
        new_query_end = orig_query_start + first_fragment_length - 1

        new_id = f"{gene_name}|{node_base}_A|{orig_ref_start}-{new_ref_end}|{orig_query_start}-{new_query_end}|{identity}|{first_fragment_length}"
        new_record = SeqRecord(Seq(first_fragment_str), id=new_id, description="")
        trimmed = trim_to_reference_bounds(new_record, reference_seq)
        if trimmed:
            sequences.append(trimmed)

    # Process last fragment
    last_fragment_str = seq_str[n_regions[-1].end():]
    if last_fragment_str:
        last_fragment_length = len(last_fragment_str)
        new_ref_start = orig_ref_end - last_fragment_length + 1
        new_query_start = orig_query_end - last_fragment_length + 1

        new_id = f"{gene_name}|{node_base}_B|{new_ref_start}-{orig_ref_end}|{new_query_start}-{orig_query_end}|{identity}|{last_fragment_length}"
        new_record = SeqRecord(Seq(last_fragment_str), id=new_id, description="")
        trimmed = trim_to_reference_bounds(new_record, reference_seq)
        if trimmed:
            sequences.append(trimmed)

    return sequences

def parse_header(record_id):
    """Parses the header format created by Step 3."""
    parts = record_id.split("|")
    gene_name = parts[0]
    node_id = parts[1]
    ref_region = tuple(map(int, parts[2].split("-")))
    query_region = tuple(map(int, parts[3].split("-")))
    identity = float(parts[4])
    alignment_length = int(parts[5])
    return gene_name, node_id, ref_region, query_region, identity, alignment_length

@lru_cache(maxsize=None)
def find_overlap(seq1, seq2, min_overlap=40):
    """Find overlaps in both directions and return best"""
    # Forward overlap (seq1 end -> seq2 start)
    fwd_overlap, fwd_start = find_overlap_single(seq1, seq2, min_overlap)
    # Reverse overlap (seq2 end -> seq1 start)
    rev_overlap, rev_start = find_overlap_single(seq2, seq1, min_overlap)

    if fwd_overlap >= rev_overlap:
        return fwd_overlap, fwd_start, "forward"
    else:
        return rev_overlap, rev_start, "reverse"

@lru_cache(maxsize=None)
def find_overlap_single(seq1, seq2, min_overlap):
    """ Find overlap between two sequences in a single direction """
    max_overlap = 0
    best_start = 0

    for i in range(min_overlap, min(len(seq1), len(seq2)) + 1):
        if seq1[-i:] == seq2[:i]:
            if i > max_overlap:
                max_overlap = i
                best_start = i

    return max_overlap, best_start # Return overlap size and start position

def can_merge_nodes(node1_id, node2_id):
    """Check if nodes can be merged, including parts from same node"""
    base1 = node1_id.split("|")[1].split("_")[0]
    base2 = node2_id.split("|")[1].split("_")[0]

    # If they are from the same NODE, check for different suffixes (A, B)
    if base1 == base2:
        suffix1 = node1_id.split("|")[1].split("_")[-1] if "_" in node1_id.split("|")[1] else ""
        suffix2 = node2_id.split("|")[1].split("_")[-1] if "_" in node2_id.split("|")[1] else ""
        return suffix1 != suffix2

    return True  # Different NODEs can always merge

def merge_with_overlaps(contigs, reference_seq, min_overlap=40, log_file=None, score_weights=(6, 4, 10)):
    """Recursive merge function to find best merge chains."""
    results = []
    merge_logs = []

    if not contigs:
        return []

    # Sort contigs by reference position to optimize merging
    trimmed_contigs = []
    for contig in contigs:
        ref_start = int(contig.id.split("|")[2].split("-")[0])
        trimmed_contigs.append((ref_start, contig))

    trimmed_contigs.sort(key=lambda x: x[0])
    trimmed_contigs = [x[1] for x in trimmed_contigs]

    merge_logs.append(f"\nAvailable contigs for merging ({len(trimmed_contigs)}):")
    for c in trimmed_contigs:
        merge_logs.append(f"  {c.id}")

    # Recursive function to chain merges
    def chain_merge(current_seq, current_nodes, used_ids, depth=0):
        if depth > 5: # Limit recursion depth to prevent infinite loops
            return

        for contig in trimmed_contigs:
            if contig.id in used_ids:
                continue

            # Check if these nodes are allowed to merge
            if not can_merge_nodes(current_nodes[-1], contig.id):
                continue

            # Try overlap in both directions
            forward_overlap = find_overlap_single(current_seq, str(contig.seq), min_overlap)
            reverse_overlap = find_overlap_single(str(contig.seq), current_seq, min_overlap)

            # Check forward merge
            if forward_overlap[0] >= min_overlap:
                new_seq = current_seq + str(contig.seq)[forward_overlap[1]:]
                new_nodes = current_nodes + [contig.id]
                new_used = used_ids | {contig.id}

                # Calculate metrics
                similarity = calculate_similarity(reference_seq, new_seq, 0, len(reference_seq))
                completeness = (len(new_seq) / len(reference_seq)) * 100
                score = calculate_score(similarity, completeness, score_weights)

                results.append({
                    'sequence': new_seq, 'nodes': new_nodes, 'similarity': similarity,
                    'completeness': completeness, 'score': score, 'type': 'merged'
                })

                # Continue chaining
                if similarity > 75: # Only continue chain if merge is good
                     chain_merge(new_seq, new_nodes, new_used, depth + 1)

            # Check reverse merge (contig comes first)
            if reverse_overlap[0] >= min_overlap:
                new_seq = str(contig.seq) + current_seq[reverse_overlap[1]:]
                new_nodes = [contig.id] + current_nodes
                new_used = used_ids | {contig.id}

                # Calculate metrics
                similarity = calculate_similarity(reference_seq, new_seq, 0, len(reference_seq))
                completeness = (len(new_seq) / len(reference_seq)) * 100
                score = calculate_score(similarity, completeness, score_weights)

                results.append({
                    'sequence': new_seq, 'nodes': new_nodes, 'similarity': similarity,
                    'completeness': completeness, 'score': score, 'type': 'merged'
                })

                # Continue chaining
                if similarity > 75:
                    chain_merge(new_seq, new_nodes, new_used, depth + 1)

    # Start a chain from each contig
    for i, start_contig in enumerate(trimmed_contigs):
        chain_merge(str(start_contig.seq), [start_contig.id], {start_contig.id})

    return results

def trim_to_reference_bounds(seq_record, reference_seq):
    """Trim sequence to match reference bounds using local alignment"""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    # Get reference ends (first/last 50bp)
    ref_start_seq = str(reference_seq[:50])
    ref_end_seq = str(reference_seq[-50:])
    seq_str = str(seq_record.seq)

    if not seq_str or not ref_start_seq or not ref_end_seq:
        return None # Cannot trim empty sequences

    # Get original coordinates from Step 3 header
    try:
        parts = seq_record.id.split("|")
        orig_ref_start, original_ref_end = map(int, parts[2].split("-"))
        orig_query_start, orig_query_end = map(int, parts[3].split("-"))
    except (IndexError, ValueError):
        print(f"Warning: Could not parse header for trimming: {seq_record.id}")
        return seq_record # Return untrimmed

    # Find start alignment
    start_idx = 0
    try:
        start_alignments = aligner.align(ref_start_seq, seq_str)
        if start_alignments:
            best_start = max(start_alignments, key=lambda a: a.score)
            if best_start.score >= 20: # Use a modest score threshold
                start_idx = best_start.aligned[1][0][0]
    except Exception as e:
        print(f"Warning: Start alignment failed for {seq_record.id}: {e}")

    # Find end alignment
    end_idx = len(seq_str)
    try:
        end_alignments = aligner.align(ref_end_seq, seq_str)
        if end_alignments:
            best_end = max(end_alignments, key=lambda a: a.score)
            if best_end.score >= 20:
                end_idx = best_end.aligned[1][-1][1]
    except Exception as e:
        print(f"Warning: End alignment failed for {seq_record.id}: {e}")

    # Trim sequence
    if start_idx < end_idx:
        trimmed_seq_str = seq_str[start_idx:end_idx]
        if not trimmed_seq_str:
            return None # Trimmed to nothing

        # Update coordinates
        new_ref_length = len(trimmed_seq_str)
        new_query_start = orig_query_start + start_idx
        new_query_end = new_query_start + new_ref_length - 1

        # Update reference start based on query trim
        new_ref_start = max(1, orig_ref_start) # Simple approximation
        new_ref_end = new_ref_start + new_ref_length - 1

        new_id = f"{parts[0]}|{parts[1]}|{new_ref_start}-{new_ref_end}|{new_query_start}-{new_query_end}|{parts[4]}|{new_ref_length}"

        return SeqRecord(Seq(trimmed_seq_str), id=new_id, description="trimmed")

    return seq_record # Return original if trimming failed

def remove_duplicate_contigs(contigs):
    """Remove duplicate contigs keeping the one with highest identity"""
    seen_sequences = {}
    duplicates_removed = 0

    for contig in contigs:
        seq_key = str(contig.seq)
        if not seq_key:
            continue

        identity = float(contig.id.split("|")[4])

        if seq_key not in seen_sequences:
            seen_sequences[seq_key] = {'contig': contig, 'identity': identity}
        else:
            if identity > seen_sequences[seq_key]['identity']:
                seen_sequences[seq_key] = {'contig': contig, 'identity': identity}
            duplicates_removed += 1

    unique_contigs = [entry['contig'] for entry in seen_sequences.values()]
    return unique_contigs, duplicates_removed

def step4_assemble_and_evaluate(input_dir: str, reference_fasta: str, output_dir: str):
    """
    Step 4: Assembles and evaluates extracted gene contigs to get the best version.
    This is the REAL function that replaces the placeholder.
    """
    print("\n🚀 Starting Step 4: Assembling and Evaluating Genes...")
    os.makedirs(output_dir, exist_ok=True)

    score_weights = (6, 4, 10) # sim_weight, comp_weight, total
    min_identity = 98 # min identity for a single contig to be considered

    try:
        references = {record.id: record.seq for record in SeqIO.parse(reference_fasta, "fasta")}
    except FileNotFoundError:
        print(f"❌ Cannot proceed with Step 4: Reference FASTA file not found at {reference_fasta}")
        return

    # Loop over each file in the input directory (e.g., "genome1_genes.fasta")
    for filename in os.listdir(input_dir):
        if filename.endswith("_genes.fasta"):
            input_fasta = os.path.join(input_dir, filename)
            genome_name = filename.replace("_genes.fasta", "")
            print(f"\n  -> Processing genome: {genome_name}")

            # --- Per-Genome Variables ---
            all_logs = [f"--- Assembly Log for {genome_name} ---"]
            best_versions = {} # {gene_name: SeqRecord}

            # 1. Group all contigs from the file by gene
            contigs_by_gene = {}
            for record in SeqIO.parse(input_fasta, "fasta"):
                try:
                    gene_name = record.id.split("|")[0]
                    if gene_name not in contigs_by_gene:
                        contigs_by_gene[gene_name] = []
                    contigs_by_gene[gene_name].append(record)
                except IndexError:
                    print(f"Warning: Skipping malformed record ID: {record.id}")

            all_logs.append(f"Found {len(contigs_by_gene)} genes to process.")

            # 2. Process each gene
            for gene_name, contigs in contigs_by_gene.items():
                if gene_name not in references:
                    all_logs.append(f"Skipping {gene_name}: Not found in reference file.")
                    continue

                all_logs.append(f"\n--- Processing gene: {gene_name} ---")
                reference_seq = references[gene_name]

                # --- Pre-process contigs: Split Ns, trim, and remove duplicates ---
                processed_contigs = []
                for contig in contigs:
                    # Split on NNNN
                    split_records = split_sequence_on_n(contig, reference_seq)
                    processed_contigs.extend(split_records)

                unique_contigs, dup_removed = remove_duplicate_contigs(processed_contigs)
                all_logs.append(f"  {len(contigs)} initial contigs -> {len(unique_contigs)} unique, trimmed fragments (Removed {dup_removed} duplicates).")

                if not unique_contigs:
                    all_logs.append("  No valid contigs remaining after pre-processing.")
                    continue

                # --- Find best single contig ---
                best_single = None
                best_single_score = -1

                for contig in unique_contigs:
                    sim = calculate_similarity(reference_seq, str(contig.seq), 0, len(reference_seq))
                    comp = (len(contig.seq) / len(reference_seq)) * 100
                    score = calculate_score(sim, comp, score_weights)

                    if score > best_single_score:
                        best_single_score = score
                        best_single = {
                            'sequence': str(contig.seq), 'similarity': sim, 'completeness': comp,
                            'score': score, 'type': 'single', 'id': contig.id
                        }

                all_logs.append(f"  Best single contig: Score={best_single_score:.2f}, Sim={best_single['similarity']:.2f}%, Comp={best_single['completeness']:.2f}%")

                # --- Attempt merging ---
                best_merge = None
                best_merge_score = -1

                if len(unique_contigs) > 1:
                    all_logs.append(f"  Attempting to merge {len(unique_contigs)} fragments...")
                    merge_results = merge_with_overlaps(unique_contigs, reference_seq, min_overlap=40, score_weights=score_weights)

                    if merge_results:
                        # Find the best merge result
                        best_merge = max(merge_results, key=lambda x: x['score'])
                        best_merge_score = best_merge['score']
                        all_logs.append(f"  Best merge result: Score={best_merge_score:.2f}, Sim={best_merge['similarity']:.2f}%, Comp={best_merge['completeness']:.2f}% (Nodes: {len(best_merge['nodes'])})")
                    else:
                        all_logs.append("  No valid merges found.")

                # --- Final Selection ---
                final_best = None

                # Prioritize a high-quality single contig
                if best_single and best_single['similarity'] >= min_identity and best_single['completeness'] >= 95:
                    final_best = best_single
                    all_logs.append(f"  Selected: High-quality single contig.")
                # Otherwise, check if a merge is better than the best single
                elif best_merge and best_merge_score > best_single_score:
                    final_best = best_merge
                    all_logs.append(f"  Selected: Merged contig (better score).")
                # Otherwise, just use the best single contig we found
                elif best_single:
                    final_best = best_single
                    all_logs.append(f"  Selected: Best single contig (merge was not better).")

                # --- Add to final list for writing ---
                if final_best:
                    gene_id = f"{gene_name}|{final_best['type']}|Score:{final_best['score']:.2f}|Sim:{final_best['similarity']:.2f}|Comp:{final_best['completeness']:.2f}"
                    nodes_str = final_best.get('id', 'N/A') if final_best['type'] == 'single' else ";".join(final_best.get('nodes', []))
                    description = f"Nodes:{nodes_str}"

                    best_versions[gene_name] = SeqRecord(
                        Seq(final_best['sequence']),
                        id=gene_id,
                        description=description
                    )

            # --- Write final FASTA file for the genome ---
            output_fasta_file = os.path.join(output_dir, f'{genome_name}_best_genes.fasta')
            final_records = list(best_versions.values())

            if final_records:
                SeqIO.write(final_records, output_fasta_file, "fasta")
                print(f"  ✅ Assembly for {genome_name} completed. Saved {len(final_records)} genes to: {os.path.basename(output_fasta_file)}")
            else:
                print(f"  ⚠️ Assembly for {genome_name} completed, but no high-quality genes were finalized.")

            # --- Write log file ---
            log_file_path = os.path.join(output_dir, f"{genome_name}_assembly_log.txt")
            with open(log_file_path, 'w') as f:
                f.write("\n".join(all_logs))
            print(f"  ℹ️ Assembly log saved to: {os.path.basename(log_file_path)}")

    print("\n✅ Step 4 completed: All genes assembled and evaluated.")


# ==============================================================================
# --- MAIN PIPELINE CONTROL ---
# ==============================================================================

def main_pipeline(args):
    """Main function to control the multi-step pipeline."""
    print("Welcome to the Gene Contig Assembly Pipeline! 1.0 🧬")

    # --- 1. Define Directory Structure based on User's --output_dir ---
    db_dir = os.path.join(args.output_dir, "1_BLAST_DB")
    blast_results_dir = os.path.join(args.output_dir, "2_BLAST_Results")
    extracted_fasta_dir = os.path.join(args.output_dir, "3_Extracted_FASTA")
    final_assembly_dir = os.path.join(args.output_dir, "4_Final_Assembly")

    # --- 2. Define the steps using the new paths ---
    steps = {
        1: ("Create BLAST Databases", lambda: step1_create_blast_databases(
                input_dir=args.scaffold_dir,
                output_dir=db_dir
             )),
        2: ("Run BLAST Search", lambda: step2_run_blast(
                db_dir=db_dir,
                query_file=args.ref_file,
                output_dir=blast_results_dir
             )),
        3: ("Extract Contig Sequences", lambda: step3_extract_contigs(
                blast_results_dir=blast_results_dir,
                scaffold_fasta_dir=args.scaffold_dir,
                output_fasta_dir=extracted_fasta_dir,
                reference_cds_file=args.ref_file
             )),
        4: ("Assemble and Evaluate Genes", lambda: step4_assemble_and_evaluate(
                input_dir=extracted_fasta_dir,
                reference_fasta=args.ref_file,
                output_dir=final_assembly_dir
             ))
    }

    # --- 3. Create Root Output Directory ---
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Using output directory: {os.path.abspath(args.output_dir)}")

    # --- 4. User Interface ---
    while True:
        print("\n--- Available Steps ---")
        for num, (name, _) in steps.items():
            print(f"[{num}] {name}")
        print("[0] Exit")

        choice = input("Enter the step number to run (e.g., 1, or 1-3) or '0' to exit: ").strip()

        if choice == '0':
            print("Exiting pipeline. Goodbye!")
            break

        try:
            if '-' in choice:
                start, end = map(int, choice.split('-'))
                steps_to_run = range(start, end + 1)
            else:
                steps_to_run = [int(choice)]

            for step_num in steps_to_run:
                if step_num in steps:
                    name, func = steps[step_num]
                    print(f"\n=============================================")
                    print(f"Starting execution of STEP {step_num}: {name}")
                    print(f"=============================================")

                    # Auto-create sub-directories as needed
                    if step_num == 1: os.makedirs(db_dir, exist_ok=True)
                    if step_num == 2: os.makedirs(blast_results_dir, exist_ok=True)
                    if step_num == 3: os.makedirs(extracted_fasta_dir, exist_ok=True)
                    if step_num == 4: os.makedirs(final_assembly_dir, exist_ok=True)

                    func() # Run the step
                else:
                    print(f"Invalid step number: {step_num}. Skipping.")

        except ValueError:
            print("Invalid input format. Please enter a number (1, 2, 3, 4) or a range (e.g., 1-3).")
        except Exception as e:
            print(f"\n❌ An unexpected error occurred: {e}")
            import traceback
            traceback.print_exc() # Print full error for debugging

# ==============================================================================
# --- SCRIPT EXECUTION ---
# ==============================================================================

if __name__ == "__main__":
    # 1. Create the parser
    parser = argparse.ArgumentParser(
        description="A multi-step pipeline to find, extract, and assemble gene contigs from scaffolds.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # 2. Add the arguments
    parser.add_argument(
        "-s", "--scaffold_dir",
        required=True,
        help="Directory containing the input scaffold FASTA files (e.g., your /Scaffolds/ folder)."
    )
    parser.add_argument(
        "-r", "--ref_file",
        required=True,
        help="FASTA file with reference gene sequences (CDS). This will be used for both the BLAST query and the assembly reference."
    )
    parser.add_argument(
        "-o", "--output_dir",
        default="pipeline_results",
        help="Main output directory. All results and sub-folders will be created here. (Default: 'pipeline_results')"
    )

    # 3. Parse the arguments
    args = parser.parse_args()

    # 4. Run the pipeline
    main_pipeline(args)

