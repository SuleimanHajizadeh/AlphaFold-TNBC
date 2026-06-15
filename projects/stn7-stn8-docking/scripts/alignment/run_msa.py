#!/usr/bin/env python3
"""
STN7 and STN8 Multiple Sequence Alignment Script
Implements a robust Star-Alignment MSA algorithm in pure Python using Biopython's PairwiseAligner.
Aligns STN7, STN8, and combined homologs and saves the results in FASTA format.
"""

import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner

def star_alignment(sequences, ref_idx=0):
    """
    Performs multiple sequence alignment using the Star Alignment algorithm.
    Uses sequences[ref_idx] as the center of the star.
    """
    if not sequences:
        return []
    
    # 1. Initialize the aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    ref_seq = sequences[ref_idx]
    ref_id = ref_seq.id
    ref_str = str(ref_seq.seq)
    
    # Align each sequence to the reference
    alignments = []
    for i, seq in enumerate(sequences):
        if i == ref_idx:
            alignments.append((ref_str, ref_str)) # Ref aligned to itself
            continue
        
        query_str = str(seq.seq)
        # Run pairwise alignment
        alns = aligner.align(ref_str, query_str)
        best_aln = alns[0] # Take the top scoring alignment
        
        # Format alignment strings
        ref_aligned, query_aligned = format_alignment_strings(best_aln, ref_str, query_str)
        alignments.append((ref_aligned, query_aligned))
    
    # 2. Merge alignments by propagating gaps in the reference sequence
    # Let's build the final aligned sequences
    aligned_lens = [len(ref_aln) for ref_aln, _ in alignments]
    
    # Track position in each alignment
    positions = [0] * len(sequences)
    
    # Final aligned sequences as lists of characters
    final_aligned = [[] for _ in range(len(sequences))]
    
    # We will step through the reference sequence.
    # At each step, if any alignment has an insertion (gap in reference),
    # we must insert a gap in the final aligned reference and all other sequences,
    # except the one that has the insertion!
    # Let's implement this step-by-step
    
    ref_pos = 0
    ref_len = len(ref_str)
    
    # State tracking: list of pointers to current index in each alignment string
    aln_ptrs = [0] * len(sequences)
    
    while ref_pos < ref_len or any(ptr < len(alignments[i][0]) for i, ptr in enumerate(aln_ptrs)):
        # Check if there is an insertion in the query (gap in reference) in any alignment
        insertion_active = False
        insertion_indices = []
        
        for i in range(len(sequences)):
            if i == ref_idx:
                continue
            
            ptr = aln_ptrs[i]
            ref_aln_str = alignments[i][0]
            
            if ptr < len(ref_aln_str) and ref_aln_str[ptr] == '-':
                # This query has an insertion relative to reference
                insertion_active = True
                insertion_indices.append(i)
        
        if insertion_active:
            # Insert characters for queries that have insertions, and gaps for others
            for i in range(len(sequences)):
                if i in insertion_indices:
                    # Append the inserted character from the query
                    char = alignments[i][1][aln_ptrs[i]]
                    final_aligned[i].append(char)
                    aln_ptrs[i] += 1
                else:
                    # Append a gap
                    final_aligned[i].append('-')
                    # For reference or non-inserting queries, do not advance pointers yet
        else:
            # Standard match/mismatch/deletion (no insertion relative to reference)
            # Advance all pointers by 1 and append characters
            for i in range(len(sequences)):
                ptr = aln_ptrs[i]
                if ptr < len(alignments[i][0]):
                    final_aligned[i].append(alignments[i][1][ptr])
                    aln_ptrs[i] += 1
                else:
                    final_aligned[i].append('-')
            
            ref_pos += 1
            
    # Convert lists back to SeqRecord objects
    aligned_records = []
    for i, seq in enumerate(sequences):
        aligned_seq_str = "".join(final_aligned[i])
        new_rec = seq.__class__(seq.seq.__class__(aligned_seq_str), id=seq.id, description=seq.description)
        aligned_records.append(new_rec)
        
    return aligned_records

def format_alignment_strings(alignment, ref_str, query_str):
    """
    Formats alignment object coordinates into aligned strings with gaps.
    """
    ref_aligned = []
    query_aligned = []
    
    # Get alignment segments
    # alignment.aligned is a tuple of two lists of coordinate ranges:
    # e.g. (((0, 3), (4, 6)), ((0, 3), (3, 5)))
    ref_coords, query_coords = alignment.aligned
    
    last_r, last_q = 0, 0
    for r_range, q_range in zip(ref_coords, query_coords):
        r_start, r_end = r_range
        q_start, q_end = q_range
        
        # Add gaps for unaligned regions in reference (insertions in query)
        gap_len_q = q_start - last_q
        if gap_len_q > 0:
            ref_aligned.append('-' * gap_len_q)
            query_aligned.append(query_str[last_q:q_start])
            
        # Add gaps for unaligned regions in query (deletions in query)
        gap_len_r = r_start - last_r
        if gap_len_r > 0:
            ref_aligned.append(ref_str[last_r:r_start])
            query_aligned.append('-' * gap_len_r)
            
        # Add aligned match/mismatch block
        ref_aligned.append(ref_str[r_start:r_end])
        query_aligned.append(query_str[q_start:q_end])
        
        last_r, last_q = r_end, q_end
        
    # Append any remaining trailing characters
    if last_r < len(ref_str):
        ref_aligned.append(ref_str[last_r:])
        query_aligned.append('-' * (len(ref_str) - last_r))
    if last_q < len(query_str):
        ref_aligned.append('-' * (len(query_str) - last_q))
        query_aligned.append(query_str[last_q:])
        
    return "".join(ref_aligned), "".join(query_aligned)

def run_alignment_for_file(input_path, output_path):
    """
    Reads FASTA, performs MSA, and writes output.
    """
    print(f"Reading input file: {input_path}")
    records = list(SeqIO.parse(input_path, "fasta"))
    if not records:
        print(f"Error: No sequences found in {input_path}")
        return
    
    print(f"Aligning {len(records)} sequences...")
    # Use first sequence (Arabidopsis) as reference center
    aligned_records = star_alignment(records, ref_idx=0)
    
    print(f"Writing aligned output to: {output_path}")
    with open(output_path, "w") as f:
        SeqIO.write(aligned_records, f, "fasta")
    print(f"Alignment length: {len(aligned_records[0].seq)}")

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fasta_dir = os.path.join(base_dir, "data/fasta")
    msa_dir = os.path.join(base_dir, "data/msa")
    os.makedirs(msa_dir, exist_ok=True)
    
    # Align STN7 homologs
    print("\n=== Aligning STN7 Homologs ===")
    run_alignment_for_file(
        os.path.join(fasta_dir, "stn7_homologs.fasta"),
        os.path.join(msa_dir, "stn7_msa.fasta")
    )
    
    # Align STN8 homologs
    print("\n=== Aligning STN8 Homologs ===")
    run_alignment_for_file(
        os.path.join(fasta_dir, "stn8_homologs.fasta"),
        os.path.join(msa_dir, "stn8_msa.fasta")
    )
    
    # Align combined all homologs
    print("\n=== Aligning All STN7 & STN8 Homologs Combined ===")
    run_alignment_for_file(
        os.path.join(fasta_dir, "all_homologs.fasta"),
        os.path.join(msa_dir, "all_msa.fasta")
    )
    
    print("\n=== Multiple Sequence Alignment Completed Successfully ===")

if __name__ == "__main__":
    main()
