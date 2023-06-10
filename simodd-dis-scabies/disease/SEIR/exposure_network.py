from math import sin, pi

from disease.general.exposure_base import Exposure


class ExposureNetworkSEIR(Exposure):
    def __init__(self, p):
        super(ExposureNetworkSEIR, self).__init__(p)
        self.phi = p['phi']
        self.iphi = 1.0 - p['phi']
        self.q = 0
        self.init_q(p)
        self.two_pi = 2 * pi / p['t_per_year']
        self.sf_amp = p['sf_amp']

    def init_q(self, p):
        self.q = p['q'] * 364.0 / p['t_per_year']

    def calc_foi_fast(self, t, ind, P, cnetwork, rng):
        """
        Calculate the probability that an individual is exposed to infection.
        """
        hh_members = P.housemates(ind)
        hh_infected = [x.ID for x in hh_members if x.state.infectious]

        # get NUMBER of infected individuals in ind's household
        # hh = self.phi * len(hh_infected)
        hh_infected_prop = len(hh_infected) / float(len(hh_members)) if hh_members else 0
        hh = self.phi * hh_infected_prop

        comm_infected = [x for x in cnetwork.G.neighbors(ind.ID) if P.I[x].state.infectious]

        # get PROPORTION of infected individuals in ind's network neighbourhood
        comm = self.iphi * (len(comm_infected) / float(cnetwork.G.degree(ind.ID)))

        combined = hh + comm

        # choose (potential) source      
        ind.hh_frac = hh / combined if combined > 0.0 else 0.0

        if combined > 0.0:
            if rng.random() < ind.hh_frac:
                ind.source = rng.sample(hh_infected, 1)[0]
                ind.hh_source = True
            else:
                ind.source = rng.sample(comm_infected, 1)[0]
                ind.hh_source = False

        # compute seasonal forcing term
        sf = (1.0 + (self.sf_amp * sin(t * self.two_pi))) if self.sf_amp > 0 else 1.0

        return self.q * sf * combined
