import gains.adaptive as adl
import pandas as pd
import salty
import gains as genetic
import unittest


class GuessIonTests(unittest.TestCase):

    def test_build_model_from_md(self):
        logs = "keras_8.0.0_salt_log.csv"
        salt_log = pd.read_csv("../../scripts/salt_logs/{}".format(logs))
        adl.build_model_from_md(salt_log, "density")
        adl.build_model_from_md(salt_log, "cpt")

    def test_calculate_minimum_distances(self):

        T = [298.1, 298.16]  # select narrow state variable ranges
        P = [101, 102]  # we will set MD simulation to 101 kPa and 298 K
        exp_data = ["cpt", "density"]
        data = salty.aggregate_data(exp_data, T=T, P=P)
        merged = salty.merge_duplicates(data)
        hull = merged.iloc[:, 2:4]
        adl.calculate_minimum_distances(hull, 1200, 600)

    def test_gaussian_pdf(self):
        # grab experimental data
        T = [297, 316]  # select narrow state variable ranges
        P = [99, 102]  # we will set MD simulation to 101 kPa and 298 K
        cpt = [[207, 3000]]
        exp_data = ["cpt"]
        cpt_data = salty.aggregate_data(exp_data, T=T, P=P, data_ranges=cpt)
        exp_data = ["density"]
        dens_data = salty.aggregate_data(exp_data, T=T, P=P)
        # calc KDEs from experimental data
        adl.gaussian_pdf(
            cpt_data.Data["Heat capacity at constant pressure, J/K/mol"])
        adl.gaussian_pdf(
            dens_data.Data["Specific density, kg/m<SUP>3</SUP>"])

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_build_model_from_md())
        genetic.Benchmark.run(self.test_calculate_minimum_distances())
        genetic.Benchmark.run(self.test_gaussian_pdf())


if __name__ == '__main__':
    unittest.main()
