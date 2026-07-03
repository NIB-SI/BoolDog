'''
Test Boolean network read/write functions.
'''
import os
import sys
import unittest

import vcr

from booldog.io import biomodels

# These tests call the live BioModels API. Real interactions are recorded
# once into the cassette files below and replayed on subsequent runs, so
# CI doesn't depend on BioModels' availability/rate limits. In CI
# (record_mode="none"), a missing/stale cassette fails loudly instead of
# silently falling through to a live request.
my_vcr = vcr.VCR(
    cassette_library_dir=os.path.join(os.path.dirname(__file__), "cassettes"),
    record_mode="none" if os.environ.get("CI") == "true" else "once",
)


class Test(unittest.TestCase):

    @my_vcr.use_cassette("test_fetch_model.yaml")
    def test_fetch_model(self):
        '''Test fetching a model from BioModels database.
        '''
        local_file_path = "./data/test_biomodels.xml"
        fp = biomodels.fetch_model(biomodels.EXAMPLE_MODEL_ID, local_file=local_file_path)
        self.assertEqual(fp, local_file_path)

    @my_vcr.use_cassette("test_fetch_model_invalid_mamo.yaml")
    def test_fetch_model_invalid_mamo(self):
        '''Test fetching a model from BioModels database with invalid modelling approach.
        '''
        with self.assertRaises(ValueError) as context:
            biomodels.fetch_model("BIOMD0000000500", check_modelling_approach=True)

        self.assertIn("Model approach", str(context.exception))


if __name__ == '__main__':
    unittest.main()
