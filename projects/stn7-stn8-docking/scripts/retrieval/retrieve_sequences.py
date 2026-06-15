#!/usr/bin/env python3
"""
STN7 and STN8 Homolog Sequence Retrieval Script
Downloads FASTA sequences from UniProt for high-confidence STN7 and STN8 homologs
across representative plant and algal species.
"""

import os
import requests

def download_fasta(accession, output_dir):
    """
    Downloads a single FASTA sequence from UniProt.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    print(f"Downloading {accession} from UniProt...")
    response = requests.get(url)
    if response.status_code == 200:
        content = response.text
        # Ensure we have a valid FASTA header
        if content.startswith(">"):
            return content
        else:
            print(f"Error: Response for {accession} does not appear to be a valid FASTA.")
            return None
    else:
        print(f"Error: Failed to fetch {accession} (HTTP {response.status_code})")
        return None

def main():
    # Define directories
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fasta_dir = os.path.join(base_dir, "data/fasta")
    os.makedirs(fasta_dir, exist_ok=True)

    # High-confidence homologs based on research phase
    stn7_homologs = {
        "Arabidopsis_thaliana": "Q9S713",
        "Oryza_sativa": "B9FLG7",
        "Zea_mays": "A0A3L6E9I0",
        "Chlamydomonas_reinhardtii": "Q84V18"  # STT7
    }

    stn8_homologs = {
        "Arabidopsis_thaliana": "Q9LZV4",
        "Oryza_sativa": "B7E5Q2",
        "Zea_mays": "A0A3L6ED31",
        "Chlamydomonas_reinhardtii": "Q84V17"  # STL1
    }

    print("=== Starting Sequence Retrieval for STN7 Homologs ===")
    stn7_sequences = []
    for species, acc in stn7_homologs.items():
        fasta_content = download_fasta(acc, fasta_dir)
        if fasta_content:
            # Clean headers to have clear descriptive labels
            lines = fasta_content.strip().split("\n")
            header = f">STN7_{species}_{acc} {lines[0][1:]}"
            sequence = "\n".join(lines[1:])
            stn7_sequences.append(f"{header}\n{sequence}\n")
    
    # Save STN7 homologs combined file
    stn7_path = os.path.join(fasta_dir, "stn7_homologs.fasta")
    with open(stn7_path, "w") as f:
        f.writelines(stn7_sequences)
    print(f"Saved STN7 homologs to: {stn7_path}")

    print("\n=== Starting Sequence Retrieval for STN8 Homologs ===")
    stn8_sequences = []
    for species, acc in stn8_homologs.items():
        fasta_content = download_fasta(acc, fasta_dir)
        if fasta_content:
            # Clean headers to have clear descriptive labels
            lines = fasta_content.strip().split("\n")
            header = f">STN8_{species}_{acc} {lines[0][1:]}"
            sequence = "\n".join(lines[1:])
            stn8_sequences.append(f"{header}\n{sequence}\n")

    # Save STN8 homologs combined file
    stn8_path = os.path.join(fasta_dir, "stn8_homologs.fasta")
    with open(stn8_path, "w") as f:
        f.writelines(stn8_sequences)
    print(f"Saved STN8 homologs to: {stn8_path}")

    # Combine all sequences for a global tree/analysis
    all_path = os.path.join(fasta_dir, "all_homologs.fasta")
    with open(all_path, "w") as f:
        f.writelines(stn7_sequences + stn8_sequences)
    print(f"\nSaved all homologs combined to: {all_path}")
    print("\n=== Sequence Retrieval Phase Completed Successfully ===")

if __name__ == "__main__":
    main()
