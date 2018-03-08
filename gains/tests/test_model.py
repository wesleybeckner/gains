from __future__ import absolute_import, division, print_function
import gains.engine as genetic
from gains.salt_generator import generate_solvent
import unittest


class GuessIonTests(unittest.TestCase):

    def test_1_model(self):
        target = 1000
        model_ID = "density_m3"
        generate_solvent(target, model_ID, heavy_atom_limit=300)

    def test_2_model(self):
        target = 1000
        model_ID = "density_m3"
        generate_solvent(target, model_ID, heavy_atom_limit=300, hits=11)

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_1_model)


if __name__ == '__main__':
    unittest.main()
