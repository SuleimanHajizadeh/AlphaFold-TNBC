import unittest
import tempfile
from pathlib import Path
from prepare_receptor import clean_pdb, write_box_config

class TestPrepareReceptor(unittest.TestCase):
    
    def setUp(self):
        # Create a temporary directory for file tests
        self.test_dir = tempfile.TemporaryDirectory()
        self.test_path = Path(self.test_dir.name)

    def tearDown(self):
        # Cleanup temporary files
        self.test_dir.cleanup()

    def test_clean_pdb_filtering(self):
        """clean_pdb should keep only kinase domain ATOMs, ignoring HETATMs and other conformers."""
        # Create a mock raw PDB file contents
        # Residue 50 is outside the kinase domain (e.g. transit peptide 1-63 for STN7)
        # Residue 100 is inside the kinase domain (64-535)
        # We also include a HETATM line (should be removed)
        # We also include a conformer B line (should be removed)
        mock_pdb_lines = [
            "ATOM      1  N   MET A  50      10.000  10.000  10.000  1.00 20.00           N\n",
            "ATOM      2  CA  MET A  50      11.000  10.000  10.000  1.00 20.00           C\n",
            "ATOM      3  N   ALA A 100      15.000  15.000  15.000  1.00 20.00           N\n",
            "ATOM      4  CA AALA A 100      16.000  15.000  15.000  1.00 20.00           C\n", # Conformer A
            "ATOM      5  CA BALA A 100      16.500  15.000  15.000  1.00 20.00           C\n", # Conformer B (should be ignored)
            "HETATM    6  O   HOH A 200      25.000  25.000  25.000  1.00 20.00           O\n", # HETATM (should be ignored)
        ]
        
        raw_pdb_file = self.test_path / "mock_raw.pdb"
        clean_pdb_file = self.test_path / "mock_clean.pdb"
        
        # Write mock raw PDB
        with open(raw_pdb_file, "w") as f:
            f.writelines(mock_pdb_lines)
            
        # Clean the PDB targeting kinase domain 64 to 535
        n_atoms = clean_pdb(raw_pdb_file, clean_pdb_file, 64, 535)
        
        # We expect only ATOM 3 and ATOM 4 to be kept:
        # ATOM 1 and 2: outside kinase domain range (residue 50 < 64)
        # ATOM 5: conformer B (ignored)
        # HETATM 6: not an ATOM line (ignored)
        self.assertEqual(n_atoms, 2)
        
        # Verify the file contents
        with open(clean_pdb_file) as f:
            lines = f.readlines()
            
        self.assertEqual(len(lines), 3) # 2 ATOM lines + END line
        self.assertTrue(lines[0].startswith("ATOM      3"))
        self.assertTrue(lines[1].startswith("ATOM      4"))
        self.assertTrue(lines[2].startswith("END"))
        
        # Check that the conformer column for A is cleaned (space instead of 'A')
        self.assertEqual(lines[1][16], ' ')

    def test_write_box_config(self):
        """write_box_config should output a correctly formatted AutoDock Vina configuration file."""
        protein_info = {
            "name": "Test Kinase",
            "organism": "Arabidopsis thaliana",
            "uniprot": "Q12345",
            "atp_center": (10, 11, 12),
            "box_center": (1.0, 2.0, 3.0),
            "box_size": (15.0, 15.0, 15.0)
        }
        
        config_file = self.test_path / "test_vina.config"
        write_box_config("TEST_PROT", protein_info, config_file)
        
        # Read contents and verify fields
        content = config_file.read_text()
        self.assertIn("center_x = 1.0", content)
        self.assertIn("center_y = 2.0", content)
        self.assertIn("center_z = 3.0", content)
        self.assertIn("size_x = 15.0", content)
        self.assertIn("size_y = 15.0", content)
        self.assertIn("size_z = 15.0", content)

if __name__ == "__main__":
    unittest.main()
