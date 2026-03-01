'''
Test Boolean network read/write functions.
'''
import sys
import unittest

from booldog.io import biomodels

class Test(unittest.TestCase):

    def test_fetch_model(self):
        '''Test fetching a model from BioModels database.
        '''
        local_file_path = "./data/test_biomodels.xml"
        fp = biomodels.fetch_model(biomodels.EXAMPLE_MODEL_ID, local_file=local_file_path)
        self.assertEqual(fp, local_file_path)

    def test_fetch_model_invalid_mamo(self):
        '''Test fetching a model from BioModels database with invalid modelling approach.
        '''
        with self.assertRaises(ValueError) as context:
            biomodels.fetch_model("BIOMD0000000500", check_modelling_approach=True)

        self.assertIn("Model approach", str(context.exception))


if __name__ == '__main__':
    unittest.main()
