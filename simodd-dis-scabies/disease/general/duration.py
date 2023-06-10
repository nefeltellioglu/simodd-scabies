import numpy as np

class DurationGeneratorNoRecCases(object):
    """
    Generates random durations for 'time in compartment' (e.g., infectious, exposed)
    according to a Gamma distribution with specified mean and shape parameter k (essentially
    equivalent to the sum of k exponential distributions)
    Addition of cases who do not receive treatment and do not recover unless MDA is provided
    """

    def __init__(self, mean, k, norecovery_cases_frac, rng):
        self.k = k
        self.theta = float(mean) / k
        self.rng = np.random.RandomState(rng.randint(0, 99999999))
        self.norecovery_cases_frac = norecovery_cases_frac
    def get_duration(self):
        if self.rng.random() < (1 - self.norecovery_cases_frac):
            return int(round(self.rng.gamma(self.k, self.theta)))
        else:
            return int(round(self.rng.gamma(self.k, self.theta * 20)))


class DurationGenerator(object):
    """
    Generates random durations for 'time in compartment' (e.g., infectious, exposed)
    according to a Gamma distribution with specified mean and shape parameter k (essentially
    equivalent to the sum of k exponential distributions)
    """

    def __init__(self, mean, k, rng):
        self.k = k
        self.theta = float(mean) / k
        self.rng = np.random.RandomState(rng.randint(0, 99999999))

    def get_duration(self):
        return int(round(self.rng.gamma(self.k, self.theta)))


class DurationGeneratorFixed(object):
    def __init__(self, mean):
        self.mean = mean

    def get_duration(self):
        return self.mean


class DurationGeneratorExpo(object):
    def __init__(self, mean, rng):
        self.mean = mean
        self.rng = np.random.RandomState(rng.randint(0, 99999999))

    def get_duration(self):
        return int(round(self.rng.exponential(self.mean)))