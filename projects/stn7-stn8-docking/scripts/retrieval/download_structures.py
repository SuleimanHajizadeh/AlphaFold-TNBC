#!/usr/bin/env python3
"""
STN7 and STN8 AlphaFold Structure Downloader Script
Queries the AlphaFold Protein Structure Database API programmatically
to retrieve and download the latest active PDB structures for STN7 and STN8.
"""

import os
import requests

def download_alphafold_structure(uniprot_id, output_dir):
    """
    Queries AlphaFold DB API to find and download the active PDB file for a UniProt ID.
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"Querying AlphaFold API for accession: {uniprot_id}...")
    
    try:
        response = requests.get(api_url)
        if response.status_code == 200:
            data = response.json()
            if data:
                prediction = data[0]
                pdb_url = prediction.get("pdbUrl")
                if pdb_url:
                    print(f"Found active PDB URL: {pdb_url}")
                    filename = f"{uniprot_id}_AlphaFold.pdb"
                    filepath = os.path.join(output_dir, filename)
                    
                    print(f"Downloading structure file...")
                    pdb_res = requests.get(pdb_url)
                    if pdb_res.status_code == 200:
                        with open(filepath, "w") as f:
                            f.write(pdb_res.text)
                        print(f"Successfully downloaded and saved: {filepath}")
                        return filepath
                    else:
                        print(f"Error: Failed to download PDB file from {pdb_url} (HTTP {pdb_res.status_code})")
                else:
                    print(f"Error: No pdbUrl found in the prediction record for {uniprot_id}")
            else:
                print(f"Error: Empty API response for {uniprot_id}")
        else:
            print(f"Error: AlphaFold API query failed (HTTP {response.status_code})")
    except Exception as e:
        print(f"Exception occurred while downloading {uniprot_id}: {str(e)}")
    
    # Fallback to direct download link using current active version v6 if API failed
    print("Attempting fallback direct download using v6 URL pattern...")
    fallback_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"
    filename = f"{uniprot_id}_AlphaFold.pdb"
    filepath = os.path.join(output_dir, filename)
    try:
        res = requests.get(fallback_url)
        if res.status_code == 200:
            with open(filepath, "w") as f:
                f.write(res.text)
            print(f"Successfully saved via fallback: {filepath}")
            return filepath
        else:
            print(f"Fallback download failed (HTTP {res.status_code})")
    except Exception as e:
        print(f"Fallback exception: {str(e)}")
    
    return None

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    structures_dir = os.path.join(base_dir, "data/structures")
    os.makedirs(structures_dir, exist_ok=True)

    targets = {
        "STN7": "Q9S713",
        "STN8": "Q9LZV4"
    }

    print("=== Starting AlphaFold Structure Download Phase ===")
    for name, uniprot_id in targets.items():
        print(f"\nProcessing {name} (UniProt: {uniprot_id})...")
        download_alphafold_structure(uniprot_id, structures_dir)
        
    print("\n=== AlphaFold Structure Download Phase Completed ===")

if __name__ == "__main__":
    main()
