from __future__ import absolute_import, division, print_function
import gains as genetic
from gains.salt_generator import generate_solvent
import unittest


class GuessIonTests(unittest.TestCase):

    def test_1_model(self):
        target = 1000
        model_ID = "density_m3"
        generate_solvent(target, model_ID)

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_1_model)


if __name__ == '__main__':
    unittest.main()
