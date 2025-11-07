from functools import lru_cache
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import itertools
from Levenshtein import distance as levenshtein_distance

@lru_cache(maxsize=None)
def calculate_similarity(reference_seq, merged_seq, ref_start, ref_end):
    """Calculate similarity using modified Levenshtein distance"""
    if not reference_seq or not merged_seq:
        return 0.0
    
    # Get the relevant part of reference sequence
    reference_region = str(reference_seq[ref_start:ref_end])
    merged_seq = str(merged_seq)
    
    # Calculate Levenshtein similarity
    max_length = max(len(reference_region), len(merged_seq))
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
    import re
    
    seq_str = str(record.seq)
    n_regions = list(re.finditer('N{5,}', seq_str))
    
    if not n_regions:
        return [trim_to_reference_bounds(record, reference_seq)]
        
    sequences = []
    gene_name = record.id.split('|')[0]
    node_base = record.id.split('|')[1]
    
    # Get original positions and identity
    orig_ref_start = int(record.id.split("|")[2].split("-")[0])
    orig_ref_end = int(record.id.split("|")[2].split("-")[1])
    orig_query_start = int(record.id.split("|")[3].split("-")[0])
    orig_query_end = int(record.id.split("|")[3].split("-")[1])
    identity = record.id.split("|")[4]
    
    # Process first fragment
    first_fragment = seq_str[:n_regions[0].start()]
    if first_fragment:
        first_fragment_length = len(first_fragment)
        # First fragment keeps original start coordinates
        new_ref_end = orig_ref_start + first_fragment_length - 1
        new_query_end = orig_query_start + first_fragment_length - 1
        
        new_id = f"{gene_name}|{node_base}_A|{orig_ref_start}-{new_ref_end}|{orig_query_start}-{new_query_end}|{identity}|{first_fragment_length}"
        new_record = SeqRecord(Seq(first_fragment), id=new_id, description="")
        sequences.append(new_record)
    
    # Process last fragment
    last_fragment = seq_str[n_regions[-1].end():]
    if last_fragment:
        last_fragment_length = len(last_fragment)
        # Last fragment keeps original end coordinates
        new_ref_start = orig_ref_end - last_fragment_length + 1
        new_query_start = orig_query_end - last_fragment_length + 1
        
        new_id = f"{gene_name}|{node_base}_B|{new_ref_start}-{orig_ref_end}|{new_query_start}-{orig_query_end}|{identity}|{last_fragment_length}"
        new_record = SeqRecord(Seq(last_fragment), id=new_id, description="")
        sequences.append(new_record)
    
    return sequences

def parse_header(record_id):
    parts = record_id.split("|")
    gene_name = parts[0]
    ref_region = tuple(map(int, parts[2].split("-")))
    query_region = tuple(map(int, parts[3].split("-")))
    identity = float(parts[4])
    alignment_length = int(parts[5])
    return gene_name, ref_region, query_region, identity, alignment_length

def identify_complete_genes(input_fasta, reference_fasta, identity_threshold=95, completeness_threshold=98):
    """Identify complete genes after trimming to reference bounds"""
    references = {record.id: record.seq for record in SeqIO.parse(reference_fasta, "fasta")}
    complete_genes = []
    problematic_genes = []
    complete_gene_names = set()

    for record in SeqIO.parse(input_fasta, "fasta"):
        gene_name = record.id.split("|")[0]
        reference_seq = references.get(gene_name)
        if not reference_seq:
            continue
            
        # Trim sequence first to reference bounds
        trimmed_record = trim_to_reference_bounds(record, reference_seq)
        
        # Check completeness with trimmed sequence and reference
        _, ref_region, query_region, identity, alignment_length = parse_header(trimmed_record.id)
        completeness = (len(trimmed_record.seq) / len(reference_seq)) * 100

        if identity >= identity_threshold and completeness >= completeness_threshold:
            complete_genes.append(trimmed_record)
            complete_gene_names.add(gene_name)
        else:
            problematic_genes.append(trimmed_record)

    problematic_genes = [record for record in problematic_genes 
                        if record.id.split("|")[0] not in complete_gene_names]
    
    return complete_genes, problematic_genes, complete_gene_names

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
    
    # Se são do mesmo NODE, verifica se têm sufixos diferentes (A, B)
    if base1 == base2:
        suffix1 = node1_id.split("|")[1].split("_")[-1]
        suffix2 = node2_id.split("|")[1].split("_")[-1]
        return suffix1 != suffix2
    
    return True  # Diferentes NODEs podem sempre mergear

def find_best_overlaps(contigs, min_overlap=40, max_overlaps=5):
    """Find best overlaps considering sequence position and overlap size for each contig"""
    overlap_graph = {}
    
    for contig1 in contigs:
        overlaps = []
        ref_start1 = int(contig1.id.split("|")[2].split("-")[0])
 
        for contig2 in contigs:
            if contig1.id != contig2.id and can_merge_nodes(contig1.id, contig2.id):
                ref_start2 = int(contig2.id.split("|")[2].split("-")[0])
                overlap, start = find_overlap_single(str(contig1.seq), str(contig2.seq), min_overlap)
                
                if overlap >= min_overlap:
                    overlaps.append({
                        'contig': contig2,
                        'overlap': overlap,
                        'start': start,
                        'ref_pos': ref_start2,
                        'score': overlap + (1 if abs(ref_start1 - ref_start2) < 100 else 0)
                    })
        
        overlap_graph[contig1.id] = sorted(overlaps, 
                                         key=lambda x: x['score'], 
                                         reverse=True)[:max_overlaps]
    
    return overlap_graph

def merge_with_overlaps(contigs, reference_seq, min_overlap=40, log_file=None, score_weights=(6, 4, 10)):
    merge_logs = []
    for contig1 in contigs:
        for contig2 in contigs:
            if contig1.id != contig2.id:
                merge_logs.append(f"\nChecking merge possibility: {contig1.id} + {contig2.id}")
                if can_merge_nodes(contig1.id, contig2.id):
                    overlap = find_overlap(str(contig1.seq), str(contig2.seq), min_overlap)
                    merge_logs.append(f"Overlap found: {overlap}")

    results = []
    merge_logs = []
    
    # First trim all contigs and sort by reference position 
    trimmed_contigs = []
    for contig in contigs:
        trimmed = trim_to_reference_bounds(contig, reference_seq)
        if trimmed:
            ref_start = int(trimmed.id.split("|")[2].split("-")[0])
            trimmed_contigs.append((ref_start, trimmed))
    
    # Sort contigs by reference position
    trimmed_contigs.sort(key=lambda x: x[0])
    trimmed_contigs = [x[1] for x in trimmed_contigs]
    
    merge_logs.append(f"\nAvailable contigs for merging (sorted by ref position):")
    for c in trimmed_contigs:
        merge_logs.append(f"  {c.id}")
    
    # Try merging adjacent contigs first
    for i in range(len(trimmed_contigs)):
        for j in range(i+1, len(trimmed_contigs)):
            contig1 = trimmed_contigs[i]
            contig2 = trimmed_contigs[j]
            
            merge_logs.append(f"\nTrying merge: {contig1.id} + {contig2.id}")
            
            # Try merge and continue chain if successful
            def chain_merge(current_seq, current_nodes, used_ids, depth=0):
                if depth > 5:
                    return
                    
                for contig in trimmed_contigs:
                    if contig.id in used_ids:
                        continue
                        
                    # Try overlap in both directions and get best
                    forward_overlap = find_overlap(current_seq, str(contig.seq), min_overlap)
                    reverse_overlap = find_overlap(str(contig.seq), current_seq, min_overlap)
                    
                    # Get best overlap
                    overlap, start, direction = max([forward_overlap, reverse_overlap], 
                                                 key=lambda x: x[0])
                    
                    if overlap >= min_overlap:
                        # Merge sequences and update nodes
                        new_seq = (current_seq + str(contig.seq)[start:] 
                                 if direction == "forward" 
                                 else str(contig.seq) + current_seq[start:])
                        new_nodes = current_nodes + [contig.id]
                        print(new_nodes)
                        
                        # Calculate metrics against reference
                        similarity = calculate_similarity(reference_seq, new_seq, 0, len(new_seq)) 
                        completeness = (len(new_seq) / len(reference_seq)) * 100
                        score = calculate_score(similarity, completeness, score_weights)
                        print(completeness, similarity, score)
                        
                        
                        # Log details for debugging
                        merge_logs.append(f"\nFound merge:")
                        merge_logs.append(f"Nodes: {' -> '.join(new_nodes)}")
                        merge_logs.append(f"Overlap: {overlap}bp")
                        merge_logs.append(f"Similarity: {similarity:.2f}%")
                        merge_logs.append(f"Completeness: {completeness:.2f}%")
                        
                        # Save all merges
                        results.append({
                            'sequence': new_seq,
                            'nodes': new_nodes,
                            'similarity': similarity,
                            'completeness': completeness,
                            'score': score,
                            'overlap': overlap
                        })
                        
                        # Continue chain if similarity > 75%
                        if similarity > 75:
                            new_used = used_ids | {contig.id}
                            chain_merge(new_seq, new_nodes, new_used, depth + 1)
            
            chain_merge(str(contig1.seq), [contig1.id], {contig1.id})
    
    # Write logs to file
    if log_file:
        with open(log_file, 'a') as f:
            f.write('\n'.join(merge_logs))
    
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
    ref_start = str(reference_seq[:50])
    ref_end = str(reference_seq[-50:])
    seq_str = str(seq_record.seq)

    # Get original coordinates
    parts = seq_record.id.split("|")
    orig_ref_start, original_ref_end = map(int, parts[2].split("-"))
    orig_query_start, orig_query_end = map(int, parts[3].split("-"))
    

    # Find start alignment
    start_alignments = aligner.align(ref_start, seq_str)
    start_idx = 0
    if start_alignments:
        best_start = start_alignments[0]
        if best_start.score >= 80:
            start_idx = best_start.aligned[1][0][0]
            
    # Find end alignment
    end_alignments = aligner.align(ref_end, seq_str)
    end_idx = len(seq_str)
    if end_alignments:
        best_end = end_alignments[0]
        if best_end.score >= 80:
            end_idx = best_end.aligned[1][-1][1]
    
    # Trim sequence
    if start_idx < end_idx:
        trimmed_seq = seq_str[start_idx:end_idx]

        # Update coordinates correctly
        new_ref_length = end_idx - start_idx

        # If trimmed at start, query_start becomes 1
        if start_idx > 0:
            new_query_start = 1
        else:
            new_query_start = orig_query_start

        # Calculate new query end based on trimming
        new_query_end = new_query_start + len(trimmed_seq) - 1

        # Update query coordinates based on trimming
        new_query_start = orig_query_start + start_idx
        new_query_end = orig_query_start + end_idx - 1
        new_ref_length = abs(end_idx - start_idx) + 1
        new_id = f"{parts[0]}|{parts[1]}|{orig_ref_start}-{orig_ref_start+new_ref_length-1}|{new_query_start}-{new_query_end}|{parts[4]}|{new_ref_length}"
        
        # Log trimming details
        #print(f"Trimming {seq_record.id}: start_idx={start_idx}, end_idx={end_idx}, trimmed_seq={trimmed_seq}")
        
        return SeqRecord(Seq(trimmed_seq), id=new_id, description="")
    
    return seq_record

def merge_problematic_genes(contigs, reference_seq, output_fasta, log_file, complete_gene_names, merge_log_file):
    """Process genes with continuous merging until completion"""
    all_merges = []
    all_logs = []

    def find_overlap_with_start_codon(seq1, seq2, min_overlap=20):
        """Special overlap finder that accounts for start codons"""
        # If first sequence starts with ATG, try matching after it
        if seq1.startswith("ATG"):
            overlap_after_start = find_overlap(seq1[3:], seq2, min_overlap)
            if overlap_after_start[0] >= min_overlap:
                return (overlap_after_start[0] + 3, overlap_after_start[1], "start_codon")
        
        # Try normal overlap
        return find_overlap(seq1, seq2, min_overlap)


    
    # First process all contigs - split Ns and trim
    processed_contigs = []
    for record in contigs:
        if 'N' in str(record.seq):
            split_records = split_sequence_on_n(record, reference_seq)
            for split_record in split_records:
                trimmed = trim_to_reference_bounds(split_record, reference_seq)
                if trimmed:
                    processed_contigs.append(trimmed)
                    all_logs.append(f"  Split and trimmed: {trimmed.id}")
        else:
            trimmed = trim_to_reference_bounds(record, reference_seq)
            if trimmed:
                processed_contigs.append(trimmed)
                all_logs.append(f"  Trimmed: {trimmed.id}")
    
    # Remove duplicates
    unique_contigs, num_removed = remove_duplicate_contigs(processed_contigs)
    all_logs.append(f"\nRemoved {num_removed} duplicate contigs")
    all_logs.append("Contigs for merging:")    
    for contig in unique_contigs:
        all_logs.append(f"  {contig.id}")

    # Now proceed with merge using cleaned contigs
    remaining_contigs = unique_contigs
    best_merge = None
    round_num = 1
    
    all_logs.append(f"\nProcessing gene with {len(unique_contigs)} contigs after splitting Ns")
    
    while remaining_contigs and round_num <= 5:
        all_logs.append(f"\nRound {round_num}:")
        merge_results = merge_with_overlaps(remaining_contigs, reference_seq, log_file=merge_log_file)
        
        if merge_results:
            current_merge = max(merge_results, key=lambda x: x['score'])
            all_logs.append(f"Best merge this round:")
            all_logs.append(f"Nodes: {' + '.join(current_merge['nodes'])}")
            all_logs.append(f"Similarity: {current_merge['similarity']:.2f}%")
            all_logs.append(f"Completeness: {current_merge['completeness']:.2f}%")
            
            if current_merge['similarity'] >= 50:
                if not best_merge or current_merge['score'] > best_merge['score']:
                    best_merge = current_merge
                    all_logs.append("New best merge found!")
                    
                    remaining_contigs = [c for c in remaining_contigs 
                                       if c.id not in current_merge['nodes']]
                    
                    if current_merge['completeness'] >= 100:
                        all_logs.append("Achieved target completeness")
                        break
        else:
            all_logs.append("No valid merges found")
            break
            
        round_num += 1
    
    # Write merge logs
    with open(merge_log_file, 'a') as f:
        f.write('\n'.join(all_logs))
    
    if best_merge:
        return [best_merge]
    return []

def identify_genes_with_criteria(input_fasta, reference_fasta, identity_threshold, completeness_threshold, exclude_genes=None):
    """Identify genes meeting specific criteria, excluding certain genes"""
    references = {record.id: record.seq for record in SeqIO.parse(reference_fasta, "fasta")}
    found_genes = []
    remaining_genes = []
    found_gene_names = set()
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        gene_name = record.id.split("|")[0]
        
        # Skip if gene should be excluded
        if exclude_genes and gene_name in exclude_genes:
            remaining_genes.append(record)
            continue
            
        ref_region = tuple(map(int, record.id.split("|")[2].split("-")))
        identity = float(record.id.split("|")[4])
        alignment_length = int(record.id.split("|")[5])
        reference_seq = references.get(gene_name)
        
        if not reference_seq:
            continue

        completeness = (alignment_length / len(reference_seq)) * 100

        if identity >= identity_threshold and completeness >= completeness_threshold:
            found_genes.append(record)
            found_gene_names.add(gene_name)
        else:
            remaining_genes.append(record)

    remaining_genes = [record for record in remaining_genes 
                      if record.id.split("|")[0] not in found_gene_names]
    
    return found_genes, remaining_genes, found_gene_names

def get_gene_quality(gene):
    """Calculate quality score for a gene based on header or description"""
    try:
        # For merged genes (check description first)
        if gene.description and "Sim:" in gene.description:
            similarity = float(gene.description.split("Sim:")[1].split("%")[0])
            completeness = float(gene.description.split("Comp:")[1].split("%")[0])
            return similarity * completeness / 100
            
        # For original genes (header format)
        parts = gene.id.split("|")
        if len(parts) >= 6:
            identity = float(parts[4])
            alignment_length = float(parts[5].split('|')[0])  # Handle possible Ex_1 suffix
            return identity * alignment_length
            
        return 0.0
        
    except (IndexError, ValueError):
        return 0.0

def filter_contigs_by_highest_identity(contigs):
    """Filter contigs to retain only the highest identity fragment for each region"""
    filtered_contigs = {}
    
    for contig in contigs:
        identity = float(contig.id.split("|")[4])
        region = contig.id.split("|")[3]  # Use the region as the key
        
        if region not in filtered_contigs or identity > filtered_contigs[region]['identity']:
            filtered_contigs[region] = {'contig': contig, 'identity': identity}
    
    return [entry['contig'] for entry in filtered_contigs.values()]

def merge_problematic_genes_relaxed(contigs, reference_seq, log_file, min_relaxed_overlap=5, score_weights=(6, 4, 10)):
    """Last resort merge attempt with relaxed criteria"""
    all_logs = []
    best_chain = None
    best_chain_score = -1
    
    # Sort contigs by reference position
    sorted_contigs = sorted(contigs, 
                          key=lambda x: int(x.id.split("|")[2].split("-")[0]))
    
    def try_chain_merge(current_seq, current_nodes, used_ids):
        nonlocal best_chain, best_chain_score
        
        # Calculate current chain metrics
        similarity = calculate_similarity(reference_seq, current_seq, 0, len(reference_seq))
        completeness = (len(current_seq) / len(reference_seq)) * 100
        
        # Add heavy penalty for over-completeness
        if completeness > 100:
            penalty = (completeness - 100) * 1  # Increased penalty for relaxed merges
            adjusted_completeness = max(0, 100 - penalty)  # Prevent negative values
        else:
            adjusted_completeness = completeness
        
        # Adjust score calculation to heavily penalize over-extension
        score = calculate_score(similarity, adjusted_completeness, score_weights)
        
        # Only consider if completeness is reasonable
        if completeness <= 200:  # Maximum 200% of reference length
            if score > best_chain_score:
                best_chain = {
                    'sequence': current_seq,
                    'nodes': current_nodes.copy(),
                    'similarity': similarity,
                    'completeness': completeness,
                    'score': score
                }
                best_chain_score = score
        
        # Try adding more contigs
        for contig in sorted_contigs:
            if contig.id in used_ids:
                continue
                
            # Very relaxed overlap criteria
            overlap = find_overlap(current_seq, str(contig.seq), min_overlap=min_relaxed_overlap)
            if overlap[0] >= min_relaxed_overlap:  # Accept smaller overlaps
                new_seq = (current_seq + str(contig.seq)[overlap[1]:] 
                          if overlap[2] == "forward"
                          else str(contig.seq) + current_seq[overlap[1]:])
                new_nodes = current_nodes + [contig.id]
                try_chain_merge(new_seq, new_nodes, used_ids | {contig.id})
    
    # Try starting chain from each contig
    for start_contig in sorted_contigs:
        try_chain_merge(str(start_contig.seq), [start_contig.id], {start_contig.id})
    
    return best_chain if best_chain and best_chain['completeness'] > 50 else None

def get_base_gene_name(gene_name):
    """Remove _extended suffix if present"""
    return gene_name.replace("_extended", "")

def remove_duplicate_contigs(contigs):
    """Remove duplicate contigs keeping the one with highest identity"""
    seen_sequences = {}
    unique_contigs = []
    duplicates_removed = 0
    
    for contig in contigs:
        # Extract full node info and coordinates from header
        parts = contig.id.split("|")
        node_full = parts[1]  # Full NODE info
        coords = f"{parts[2]}|{parts[3]}"  # Both reference and query coordinates
        seq_key = (str(contig.seq), node_full, coords)
        
        # Get identity from header
        identity = float(parts[4])
        
        if seq_key not in seen_sequences:
            seen_sequences[seq_key] = {'contig': contig, 'identity': identity}
        else:
            # If we find a duplicate with higher identity, replace it
            if identity > seen_sequences[seq_key]['identity']:
                seen_sequences[seq_key] = {'contig': contig, 'identity': identity}
            duplicates_removed += 1
    
    # Convert dictionary values back to list
    unique_contigs = [entry['contig'] for entry in seen_sequences.values()]
            
    return unique_contigs, duplicates_removed

def process_files_in_directory(input_dir, reference_fasta, output_dir, 
                             score_weights=(6, 4, 10),
                             min_overlap=40,
                             min_relaxed_overlap=5,
                             min_dup_length=100,min_identity=98):
    """Process all files in the input directory"""
    
    # Load references
    references = {record.id: record.seq for record in SeqIO.parse(reference_fasta, "fasta")}
    
    for filename in os.listdir(input_dir):
        if filename.endswith("_genes.fasta"):
            input_fasta = os.path.join(input_dir, filename)
            genome_name = filename.replace("_scaffolds_genes.fasta", "")
            print(f"Processing: {genome_name}")
            
            # Initialize for each genome
            best_versions = {}
            passed_genes = []
            failed_genes = []
            all_logs = []
            details_logs = []
            merged_genes = []
            found_gene_names = set()
            
            # Group contigs by gene
            contigs_by_gene = {}
            for record in SeqIO.parse(input_fasta, "fasta"):
                gene_name = record.id.split("|")[0]
                if gene_name not in contigs_by_gene:
                    contigs_by_gene[gene_name] = []
                contigs_by_gene[gene_name].append(record)
            
            # Process each gene
            for gene_name, contigs in contigs_by_gene.items():
                print(f"Processing gene: {gene_name}")
                all_logs.append(f"\nProcessing gene: {gene_name}")
                all_logs.append("="*50)  # Add separator
                
                # Filter contigs by minimum identity first
                filtered_contigs = []
                all_logs.append("\nInitial contigs:")
                for contig in contigs:
                    identity = float(contig.id.split("|")[4])
                    all_logs.append(f"  {contig.id}")
                    if identity >= min_identity:
                        filtered_contigs.append(contig)
                    else:
                        all_logs.append(f"    Discarded: Identity {identity:.2f}% below threshold {min_identity}%")
                
                all_logs.append(f"\nRetained {len(filtered_contigs)} contigs after identity filtering")
                
                if not filtered_contigs:
                    all_logs.append("No contigs passed identity threshold, skipping gene")
                    continue


                # Continue processing with filtered_contigs instead of contigs
                copies, has_duplicates = detect_and_process_duplicates(
                    filtered_contigs,  # Use filtered contigs
                    references[gene_name],
                    min_identity=min_identity,  # Use same threshold
                    min_length=min_dup_length
                )

                if has_duplicates:
                    # Log duplicates but don't add to best_versions yet
                    all_logs.append(f"\nDetected potential duplicates for {gene_name}")
                    details_logs.append(f"\nDetected potential duplicates for {gene_name}")
                    for copy in copies:
                        details_logs.append("\nEvidence:")
                        details_logs.append(f"Nodes: {' ; '.join(copy['nodes'])}")
                        details_logs.append(f"Similarity between nodes: {copy['similarity']:.2f}%")
                        details_logs.append(f"Corresponding region size: {copy['overlap_size']} bp")
                        details_logs.append(f"Minimum number of possible copies: {copy['num_copies']}")
                        details_logs.append("Overlapping sequences:")
                        details_logs.append(f"Positional info: {copy['positional_info']}")
                        details_logs.append(f"Identity: {copy['identity1']:.2f}% and {copy['identity2']:.2f}%")
                        details_logs.append(f"Region1: {copy['overlap_regions']['region1']}")
                        details_logs.append(f"Region2: {copy['overlap_regions']['region2']}")

                # Process gene normally regardless of duplicates
                # Try to find complete gene first
                best_single = None
                best_single_score = -1
                
                # Process only filtered contigs from here on
                all_logs.append("\nTrying single contigs:")
                for contig in filtered_contigs:
                    trimmed = trim_to_reference_bounds(contig, references[gene_name])
                    if trimmed:
                        sim = calculate_similarity(references[gene_name], str(trimmed.seq), 0, len(trimmed.seq))
                        comp = (len(trimmed.seq) / len(references[gene_name])) * 100
                        score = calculate_score(sim, comp, score_weights)
                        
                        all_logs.append(f"\nEvaluating: {trimmed.id}")
                        all_logs.append(f"  Similarity: {sim:.2f}%")
                        all_logs.append(f"  Completeness: {comp:.2f}%")
                        all_logs.append(f"  Score: {score:.2f}")
                        
                        if score > best_single_score:
                            best_single_score = score
                            best_single = {
                                'sequence': str(trimmed.seq),
                                'similarity': sim,
                                'completeness': comp,
                                'score': score,
                                'type': 'single',
                                'id': trimmed.id
                            }
                            all_logs.append(" New best single contig!")

                if best_single:
                    all_logs.append(f"\nBest single contig found:")
                    all_logs.append(f"  ID: {best_single['id']}")
                    all_logs.append(f"  Score: {best_single_score:.2f}")
                    all_logs.append(f"  Similarity: {best_single['similarity']:.2f}%")
                    all_logs.append(f"  Completeness: {best_single['completeness']:.2f}%")

                # Try merging if needed
                best_merge = None
                if not (best_single and best_single['completeness'] >= 95 and best_single['similarity'] >= 95):
                    all_logs.append("\nSingle contig doesn't meet quality criteria. Attempting merge strategies...")

                    # First try normal merge
                    all_logs.append("\nTrying standard merge...")
                    merged_result = merge_problematic_genes(
                        filtered_contigs,
                        references[gene_name],
                        os.path.join(output_dir, f"{genome_name}_merged.fasta"),
                        os.path.join(output_dir, f"{genome_name}.log"),
                        found_gene_names,
                        os.path.join(output_dir, f"{genome_name}_merge_details.log")
                    )
                    
                    if merged_result:
                        best_merge = max(merged_result, key=lambda x: x['score'])
                        best_merge['type'] = 'merged'
                        all_logs.append(f"\nSuccessful merge found:")
                        all_logs.append(f"  Nodes: {' -> '.join(best_merge['nodes'])}")
                        all_logs.append(f"  Similarity: {best_merge['similarity']:.2f}%")
                        all_logs.append(f"  Completeness: {best_merge['completeness']:.2f}%")
                        all_logs.append(f"  Score: {best_merge['score']:.2f}")
                    # If still not complete/correct, try relaxed merge
                    if not (best_merge and best_merge['completeness'] >= 95 and best_merge['similarity'] >= 90):
                        all_logs.append("\nStandard merge not sufficient. Trying relaxed merge...")
                        relaxed_result = merge_problematic_genes_relaxed(
                            contigs,
                            references[gene_name],
                            os.path.join(output_dir, f"{genome_name}_merge_details.log"),
                            min_relaxed_overlap=min_relaxed_overlap
                        )
                        
                        if relaxed_result:
                            if not best_merge or relaxed_result['score'] > best_merge['score']:
                                best_merge = relaxed_result
                                best_merge['type'] = 'merged_relaxed'
                                all_logs.append(f"\nSuccessful relaxed merge found:")
                                all_logs.append(f"  Nodes: {' -> '.join(best_merge['nodes'])}")
                                all_logs.append(f"  Similarity: {best_merge['similarity']:.2f}%")
                                all_logs.append(f"  Completeness: {best_merge['completeness']:.2f}%")
                                all_logs.append(f"  Score: {best_merge['score']:.2f}")   

                        # Add to merged_genes
                        seq_record = SeqRecord(
                            Seq(best_merge['sequence']),
                            id=f"{gene_name}|merged",
                            description=f"Sim:{best_merge['similarity']:.2f}%_Comp:{best_merge['completeness']:.2f}%"
                        )
                        merged_genes.append(seq_record)

                # Select best version (single or merged)
                best_version = None
                if best_merge and best_single:
                    best_version = best_merge if best_merge['score'] > best_single_score else best_single
                elif best_merge:
                    best_version = best_merge
                elif best_single:
                    best_version = best_single
                
                # Final decision logging
                all_logs.append("\nFinal decision:")
                if best_merge and best_single:
                    if best_merge['score'] > best_single_score:
                        all_logs.append("Selected: Merged version (better score)")
                    else:
                        all_logs.append("Selected: Single contig version (better score)")
                elif best_merge:
                    all_logs.append("Selected: Merged version (no good single contig)")
                elif best_single:
                    all_logs.append("Selected: Single contig version (no successful merge)")
                else:
                    all_logs.append("No valid solution found")

                all_logs.append("\n" + "="*50 + "\n")  # Add separator between genes

                # Add best version to best_versions
                if best_version:
                    seq_record = SeqRecord(
                        Seq(best_version['sequence']),
                        id=f"{gene_name}|{best_version['type']}",
                        description=f"Sim:{best_version['similarity']:.2f}%_Comp:{best_version['completeness']:.2f}%"
                    )
                    best_versions[gene_name] = seq_record

            # Write all output files
            output_fasta = os.path.join(output_dir, f"{genome_name}_genes_final.fasta")
            merge_log_file = os.path.join(output_dir, f"{genome_name}_merge_details.log")
            details_file = os.path.join(output_dir, f"{genome_name}_details.txt")
            merged_genes_file = os.path.join(output_dir, f"{genome_name}_merged_genes.fasta")
            
            # Write files
            with open(output_fasta, "w") as f:
                SeqIO.write(best_versions.values(), f, "fasta")
            
            with open(merge_log_file, "w") as f:
                f.write('\n'.join(all_logs))
                
            with open(details_file, "w") as f:
                f.write('\n'.join(details_logs))
                
            with open(merged_genes_file, "w") as f:
                SeqIO.write(merged_genes, f, "fasta")
            
            # Clear cache
            calculate_similarity.cache_clear()
            find_overlap.cache_clear()
            find_overlap_single.cache_clear()

def detect_and_process_duplicates(contigs, reference_seq, min_identity=98.00, min_length=100):
    """
    Detect and process potential gene duplicates.
    Criteria:
    - Identity to reference >= 98%
    - Sequence length >= 100bp
    - Overlapping region >= 40bp
    - Similarity between contigs < 100% in overlapping region but > 95%
    """

    #1. First filter contigs by minimum identity STRICTLY
    filtered_contigs = []
    for contig in contigs:
        identity = float(contig.id.split("|")[4])
        if identity >= min_identity:  
            filtered_contigs.append(contig)
            print(f"  {contig.id} passed identity filter")

    if len(filtered_contigs) < 2:  # Se não tivermos pelo menos 2 contigs com alta identidade
        return [], False   

    # 2. Split on N regions
    split_contigs = []
    for record in filtered_contigs:
        if 'N' in str(record.seq):
            split_records = split_sequence_on_n(record, reference_seq)
            split_contigs.extend(split_records)
        else:
            split_contigs.append(record)
    
    # 3. Trim to reference bounds
    trimmed_contigs = []
    for record in split_contigs:
        trimmed = trim_to_reference_bounds(record, reference_seq)
        if trimmed:
            trimmed_contigs.append(trimmed)
    
    # 4. Remove duplicates
    unique_contigs, _ = remove_duplicate_contigs(trimmed_contigs)
    sorted_contigs = sorted(unique_contigs, 
                          key=lambda x: int(x.id.split("|")[2].split("-")[0]))   
    
    # 5. Now detect duplicates between remaining unique contigs
    duplicate_groups = []
    for i, contig1 in enumerate(sorted_contigs):
        ref_start1, ref_end1 = map(int, contig1.id.split("|")[2].split("-"))
        
        for contig2 in sorted_contigs[i+1:]:
            ref_start2, ref_end2 = map(int, contig2.id.split("|")[2].split("-"))
            
            # Check if regions overlap in reference
            overlap_start = max(ref_start1, ref_start2)
            overlap_end = min(ref_end1, ref_end2)
            overlap_length = overlap_end - overlap_start
            
            if overlap_length >= min_length:
                identity1 = float(contig1.id.split("|")[4])
                identity2 = float(contig2.id.split("|")[4])
                
                #if identity1 >= min_identity and identity2 >= min_identity:
                    # Extract overlapping regions
                pos1_start = overlap_start - ref_start1
                pos1_end = pos1_start + overlap_length
                pos2_start = overlap_start - ref_start2
                pos2_end = pos2_start + overlap_length
                    
                overlap_region1 = str(contig1.seq)[pos1_start:pos1_end]
                overlap_region2 = str(contig2.seq)[pos2_start:pos2_end]
                    
                # Calculate similarity
                similarity = calculate_similarity(overlap_region1, overlap_region2, 0, len(overlap_region1))
                    
                if 95 <= similarity < 100:
                    duplicate_info = {
                        'nodes': [contig1.id, contig2.id],
                        'similarity': similarity,
                        'overlap_length': overlap_length,
                        'num_copies': 2,  # Sempre 2 neste caso
                        'completeness': (len(str(contig1.seq)) / len(reference_seq)) * 100,
                        'sequence': str(contig1.seq),
                        'positional_info': [pos1_start+1, pos1_end, pos2_start+1, pos2_end],
                        'identity1': identity1,
                        'identity2': identity2,
                        'overlap_regions': {
                            'region1': overlap_region1,
                            'region2': overlap_region2,
                        }
                    }
                    duplicate_groups.append(duplicate_info)
    
    if not duplicate_groups:
        return [], False
        
    # Format results
    evidence = []
    for group in duplicate_groups:
        evidence.append({
            'nodes': group['nodes'],
            'similarity': group['similarity'],
            'overlap_size': group['overlap_length'],
            'num_copies': len(group['nodes']),
            'completeness': group['completeness'],
            'sequence': group['sequence'],
            'positional_info': group['positional_info'],
            'identity1': group['identity1'],
            'identity2': group['identity2'],
            'overlap_regions': group['overlap_regions'],
        })
    
    return evidence, bool(evidence)

def merge_duplicated_gene(contigs, reference_seq, used_contigs, copy_name):
    """Special merge function for duplicated genes that allows reuse of contigs"""
    best_chain = None
    best_score = -1
    
    for start_contig in contigs:
        current_seq = str(start_contig.seq)
        current_nodes = [start_contig.id]
        
        chain_result = chain_merge_duplicates(
            current_seq,
            current_nodes,
            contigs,
            reference_seq,
            used_contigs
        )
        
        if chain_result and chain_result['score'] > best_score:
            best_chain = chain_result
            best_score = chain_result['score']
    
    if best_chain:
        # Update the chain's ID to include copy name
        best_chain['id'] = f"{best_chain['id']}_{copy_name}|merged"
    
    return best_chain

def chain_merge_duplicates(current_seq, current_nodes, contigs, reference_seq, used_contigs, 
                         min_overlap=20, depth=0, score_weights=(6, 4, 10)):
    """Chain merge for duplicated genes with relaxed constraints"""
    if depth > 10:  # Prevent infinite recursion
        return None
    
    # Get gene name from first node
    gene_name = current_nodes[0].split("|")[0]
    
    best_result = {
        'sequence': current_seq,
        'nodes': current_nodes,
        'similarity': calculate_similarity(reference_seq, current_seq, 0, len(current_seq)),
        'completeness': (len(current_seq) / len(reference_seq)) * 100,
        'id': gene_name  # Add ID to the result dictionary
    }
    best_result['score'] = calculate_score(best_result['similarity'], best_result['completeness'], score_weights)
    
    # Rest of the function remains the same...
    return best_result

def generate_summary_report(output_dir, output_summary="gene_status_summary.csv"):
    """Generate summary report in matrix format with genes as columns"""
    all_results = {}
    all_genes = set()
    
    # First pass: collect results from final fasta files
    for filename in os.listdir(output_dir):
        if filename.endswith("_genes_final.fasta"):
            genome_name = filename.replace("_genes_final.fasta", "")
            fasta_path = os.path.join(output_dir, filename)
            
            # Initialize results for this genome
            all_results[genome_name] = {
                'passed': set(),
                'failed': set(),
                'duplicated': set()
            }
            
            # Get genes from final fasta
            for record in SeqIO.parse(fasta_path, "fasta"):
                gene_name = record.id.split("|")[0]
                sim = float(record.description.split("Sim:")[1].split("%")[0])
                comp = float(record.description.split("Comp:")[1].split("%")[0])
                
                all_genes.add(gene_name)
                if sim >= 90 and comp >= 95:
                    all_results[genome_name]['passed'].add(gene_name)
                else:
                    all_results[genome_name]['failed'].add(gene_name)
            
            # Check for duplicates in details file
            details_file = os.path.join(output_dir, f"{genome_name}_details.txt")
            if os.path.exists(details_file):
                with open(details_file, 'r') as f:
                    content = f.read()
                    for line in content.split("\n"):
                        if "Detected potential duplicates for" in line:
                            gene = line.split("for ")[1].strip()
                            all_results[genome_name]['duplicated'].add(gene)
    
    # Sort genes alphabetically
    ordered_genes = sorted(list(all_genes))
    
    # Write report
    with open(os.path.join(output_dir, output_summary), 'w') as f:
        # Write header
        f.write("Genome\t" + "\t".join(ordered_genes) + "\n")
        
        # Write data for each genome
        for genome in sorted(all_results.keys()):
            row = [genome]
            for gene in ordered_genes:
                if gene in all_results[genome]['passed']:
                    status = "Complete (Duplicated)" if gene in all_results[genome]['duplicated'] else "Complete"
                elif gene in all_results[genome]['failed']:
                    status = "Failed (Duplicated)" if gene in all_results[genome]['duplicated'] else "Failed"
                else:
                    status = "Failed"  # Default for genes not found
                row.append(status)
            
            f.write("\t".join(row) + "\n")
    
    print(f"\nSummary report generated: {os.path.join(output_dir, output_summary)}")

# Add to main script
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="""
Process gene sequences with configurable parameters.

Examples:
    # Basic usage with required parameters:
    python Merge_user_friendly.py -i input_folder -r reference.fasta -o output_folder

    # Advanced usage with custom parameters:
    python Merge_user_friendly.py -i input_folder -r reference.fasta -o output_folder -s 8 -c 2 -m 30 -M 10 -l 80 -I 95

Note: total_weight is automatically calculated as sim_weight + comp_weight
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-i', '--input_dir', required=True,
                      help='Directory containing input FASTA files')
    parser.add_argument('-r', '--reference_fasta', required=True,
                      help='Path to reference FASTA file')
    parser.add_argument('-o', '--output_dir', required=True,
                      help='Directory for output files')
    
    # Optional arguments with defaults
    parser.add_argument('-s', '--sim_weight', type=int, default=6,
                      help='Weight for similarity in score calculation (default: 6)')
    parser.add_argument('-c', '--comp_weight', type=int, default=4,
                      help='Weight for completeness in score calculation (default: 4)')
    parser.add_argument('-m', '--min_overlap', type=int, default=40,
                      help='Minimum overlap for main merge strategy (default: 40)')
    parser.add_argument('-M', '--min_relaxed_overlap', type=int, default=5,
                      help='Minimum overlap for relaxed merge strategy (default: 5)')
    parser.add_argument('-l', '--min_length', type=int, default=100,
                      help='Minimum alignment length for duplicate detection (default: 100)')
    parser.add_argument('-I', '--min_identity', type=float, default=98,
                      help='Minimum identity threshold for contigs (default: 98)')
    args = parser.parse_args()

    # Process files with user parameters
    process_files_in_directory(
        input_dir=args.input_dir,
        reference_fasta=args.reference_fasta,
        output_dir=args.output_dir,
        score_weights=(args.sim_weight, args.comp_weight, args.sim_weight + args.comp_weight),
        min_overlap=args.min_overlap,
        min_relaxed_overlap=args.min_relaxed_overlap,
        min_dup_length=args.min_length,
        min_identity=args.min_identity 
    )
    
    # Generate summary report
    generate_summary_report(args.output_dir)

