from math import sin, pi

from disease.general.exposure_base import Exposure


class ExposureSEIRadd(Exposure):
    def __init__(self, p):
        """
        Initialise a new Exposure object -- calculates FOI

        :param p: parameter dictionary
        """
        super(ExposureSEIRadd, self).__init__(p)

        # transmission coefficients
        self.q = p['q'] * 364.0 / p['t_per_year']
        self.q_h = p['q_h'] * 364.0 / p['t_per_year']

        # household transmission:
        # alpha = 0 : density dependent household mixing
        # alpha = 1 : frequency dependent household mixing
        self.alpha = p.get('alpha', 1.0)

        # relative susceptibility of adults
        self.rel_susc_adult = p.get('rel_susc_adult', 1.0)

        # seasonal forcing: amplitude
        self.sf_amp = p.get('sf_amp', 0.0)
        # seasonal forcing: pre-compute 2*pi
        self.two_pi = 2 * pi / p['t_per_year']

    def calc_foi_fast(self, t, ind, pop, comm_precomp, rng):
        """
        Calculates force of infection acting on an individual based on prevalence
        of infection in household and community

        :param t: current timestep
        :param ind: current individual
        :param pop: population
        :param comm_precomp: pre-computed age-specific community FOI
        :return: FOI acting on ind
        """

        # list of hh_members
        hh_members = pop.groups['household'][ind.groups['household']]

        # list of infectious hh members (I_H)
        hh_i = sum(1 for x in hh_members if x.state.infectious)
        hh_infected = [x.ID for x in hh_members if x.state.infectious]
        N_sub_1 = len(hh_members) - 1

        # proportion/count of household members infectious (I_H / (N_H - 1)^alpha)
        hh_I_prop = hh_i / pow(N_sub_1, self.alpha) if N_sub_1 else 0

        # calculate household contribution
        hh = self.q_h * hh_I_prop

        # calculate community contribution
        # (uses pre-computed values for \sum_j (\eta_{ij} I_j / N_j) )
        comm = self.q * comm_precomp['I'][ind.age]

        # combined exposure is sum of household and community
        combined = hh + comm

        # apply relative adult susceptibility
        if ind.age > 18:
            combined *= self.rel_susc_adult

        # determine if household source and if
        #ind.hh_frac = hh / (hh + comm) if (hh + comm) > 0.0 else 0.0
        #display(hh)
        #display(hh + comm)
        ind.hh_frac = hh / combined if combined > 0.0 else 0.0
        #display(ind.hh_frac)
        #if hh > 0: print (hh / (hh+comm)), (hh / combined)

        if ind.hh_frac > 0.0 and rng.random() < ind.hh_frac:
            ind.hh_source = 1
            ind.source = rng.sample(list(hh_infected), 1)[0]
            #print ind.hh_frac, ind.hh_source, ind.source
        else:
            ind.hh_source = 0
            ind.source = 0

        # compute seasonal forcing term
        sf = (1.0 + (self.sf_amp * sin(t * self.two_pi))) if self.sf_amp > 0.0 else 1.0

        return sf * combined
