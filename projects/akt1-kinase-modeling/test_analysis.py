import unittest
import tempfile
import json
import os
import numpy as np
from pathlib import Path
from analysis import find_rank001_json, load_plddt, colour_by_band

class TestAlphaFoldAnalysis(unittest.TestCase):

    def setUp(self):
        # Create a temporary directory for file tests
        self.test_dir = tempfile.TemporaryDirectory()
        self.test_path = Path(self.test_dir.name)

    def tearDown(self):
        # Cleanup temporary files
        self.test_dir.cleanup()

    def test_find_rank001_json_not_found(self):
        """find_rank001_json should raise FileNotFoundError when no matching file is found."""
        with self.assertRaises(FileNotFoundError):
            find_rank001_json(str(self.test_path))

    def test_find_rank001_json_success(self):
        """find_rank001_json should correctly locate files matching the pattern."""
        dummy_file = self.test_path / "model_scores_rank_001_prediction.json"
        dummy_file.write_text("{}")
        
        found = find_rank001_json(str(self.test_path))
        self.assertEqual(os.path.abspath(found), os.path.abspath(dummy_file))

    def test_load_plddt(self):
        """load_plddt should correctly parse pLDDT and PTM fields from the JSON file."""
        mock_data = {
            "plddt": [92.5, 85.0, 60.0, 45.2],
            "ptm": 0.825
        }
        
        json_file = self.test_path / "mock_scores.json"
        with open(json_file, "w") as f:
            json.dump(mock_data, f)
            
        plddt, ptm = load_plddt(str(json_file))
        
        np.testing.assert_array_equal(plddt, np.array([92.5, 85.0, 60.0, 45.2]))
        self.assertEqual(ptm, 0.825)

    def test_colour_by_band(self):
        """colour_by_band should map score lists to correct standard EBI confidence color bands."""
        # Very high (>=90) -> #0053D6
        # Confident (70–90) -> #65CBF3
        # Low (50–70) -> #FFDB13
        # Very low (<50) -> #FF7D45
        scores = [95.0, 90.0, 89.9, 70.0, 69.9, 50.0, 49.9, 10.0]
        expected_colors = [
            "#0053D6", "#0053D6", 
            "#65CBF3", "#65CBF3", 
            "#FFDB13", "#FFDB13", 
            "#FF7D45", "#FF7D45"
        ]
        
        mapped_colors = colour_by_band(scores)
        self.assertEqual(mapped_colors, expected_colors)

if __name__ == "__main__":
    unittest.main()
