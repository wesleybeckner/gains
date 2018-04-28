from __future__ import absolute_import, division, print_function
from gains.salt_generator import generate_solvent
import unittest


class GuessIonTests(unittest.TestCase):

    def test_1_model(self):
        target = 1300
        model_ID = "density"
        generate_solvent(target, model_ID, heavy_atom_limit=300)


if __name__ == '__main__':
    unittest.main()
