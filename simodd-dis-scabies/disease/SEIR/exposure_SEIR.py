from math import sin, pi

from disease.general.exposure_base import Exposure


class ExposureSEIR(Exposure):
    def __init__(self, p):
        super(ExposureSEIR, self).__init__(p)
        self.phi = p['phi']
        self.iphi = 1.0 - p['phi']
        self.q = 0
        self.init_q(p)
        self.two_pi = 2 * pi / p['t_per_year']
        self.sf_amp = p['sf_amp']
        # alpha = 0 ---> density dependent household mixing
        # alpha = 1 ---> frequency dependent household mixing
        # default to density dependant mixing if alpha not specified
        self.alpha = p['alpha'] if 'alpha' in p else 0

    def init_q(self, p):
        self.q = p['q'] * 364.0 / p['t_per_year']

    def calc_foi_fast(self, t, ind, P, comm_precomp, rng):
        """
        Calculate the probability that an individual is exposed to infection.

        Currently:
            - household mixing is density dependent
            - community mixing is frequency dependent

        Options:
        1. either implement frequency dependent household mixing separately
        (if more efficient that way), or use alpha to parameterise. (**DONE**)

        2. enable sex-dependent susceptibility (possibly in a new class,
        via a check and weighting parameter.

        3. sex-dependent infectiousness would seem most easily handled via
        an alternative disease model, in which entry into I_m or I_f is
        determined by sex, enabling re-use of current 'disease_base.by_age',
        as per John's SEIAR model.
        """
        # calculate weighted household contribution to effective contact rate
        # NB: hh_members DOES include ind (so need to subtract 1 below)
        hh_members = P.groups['household'][ind.groups['household']]
        hh_infected = sum(1 for x in hh_members if x.state.infectious)
        # hh_infected = [x.ID for x in hh_members if x.state.infectious]

        # get NUMBER of infected individuals in ind's household
        #hh = self.phi * len(hh_infected)
        N_sub_1 = len(hh_members) - 1

        # get the PROPORTION OF infected individuals in ind's household
        hh_infected_prop = hh_infected / pow(N_sub_1, self.alpha) if N_sub_1 else 0.0
        hh = self.phi * hh_infected_prop

        # calculate weighted community contribution from precomputed values
        comm = self.iphi * comm_precomp['I'][ind.age]

        # compute combined effective contact rate
        combined = hh + comm

        ind.hh_frac = hh / combined if combined > 0.0 else 0.0

        ind.hh_source = 0
        ind.source = 0

        # if combined > 0.0:
        #     ind.hh_source = rng.random() < ind.hh_frac
        #
        # if ind.hh_source:
        #     ind.source = rng.sample(hh_infected, 1)[0]

        # compute seasonal forcing term
        sf = (1.0 + (self.sf_amp * sin(t * self.two_pi))) if self.sf_amp > 0 else 1.0

        return self.q * sf * combined
