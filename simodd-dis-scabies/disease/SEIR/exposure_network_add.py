from math import sin, pi

from disease.general.exposure_base import Exposure


class ExposureNetworkSEIRadd(Exposure):
    def __init__(self, p):
        super(ExposureNetworkSEIRadd, self).__init__(p)
        self.q = 0
        self.q_h = 0
        self.init_q(p)
        # for seasonal forcing
        self.two_pi = 2 * pi / p['t_per_year']
        self.sf_amp = p['sf_amp']
        # alpha = 0 ---> density dependent household mixing
        # alpha = 1 ---> frequency dependent household mixing
        # default to density dependant mixing if alpha not specified
        self.alpha = p['alpha'] if 'alpha' in p else 0

    def init_q(self, p):
        self.q = p['q'] * 364.0 / p['t_per_year']
        self.q_h = p['q_h'] * 364.0 / p['t_per_year']

    def calc_foi_fast(self, t, ind, P, cnetwork, rng):
        """
        Calculate the probability that an individual is exposed to infection.
        """
        hh_members = P.housemates(ind)
        hh_infected = [x.ID for x in hh_members if x.state.infectious]

        # get the PROPORTION OF infected individuals in ind's household
        hh_infected_prop = len(hh_infected) / pow(float(len(hh_members)), self.alpha) if hh_members else 0
        hh = self.q_h * hh_infected_prop

        comm_infected = [x for x in cnetwork.G.neighbors(ind.ID) if P.I[x].state.infectious]

        # get PROPORTION of infected individuals in ind's network neighbourhood
        comm = self.q * (len(comm_infected) / float(cnetwork.G.degree(ind.ID)))

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
