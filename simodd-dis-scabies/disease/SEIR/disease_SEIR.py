from disease.general.disease_base import DiseaseBase
from disease.general.duration import DurationGenerator, DurationGeneratorNoRecCases
from disease.general.contact_network import ContactNetwork
from disease.SEIR.exposure_SEIR_add import ExposureSEIRadd
from disease.SEIR.exposure_network import ExposureNetworkSEIR
from disease.SEIR.exposure_network_add import ExposureNetworkSEIRadd
from disease.SEIR.states_SEIR import *


# Class to do either SEIR or SIR outbreak

class DiseaseSEIR(DiseaseBase):
    def __init__(self, p, cmatrix, rng, fname, mode):
        super(DiseaseSEIR, self).__init__(p, None, cmatrix, fname, mode, basic_infection='I', disease_states=('E','I'))
        self.exposure = ExposureSEIRadd(p)
        rem = Removed(3)
        inf = Infected(2, DurationGenerator(
            p['infective_mean'] / (364.0 / p['t_per_year']), p['k'], rng))

        if p['exposed_mean'] > 0:
            exp_ = Exposed(1, DurationGenerator(
                p['exposed_mean'] / (364.0 / p['t_per_year']), p['k'], rng))
            sus = Susceptible(0)
            self.add_states(sus, exp_, inf, rem)
            self.disease_states = ['I', 'E']
        else:
            sus = Susceptible(0)
            self.add_states(sus, inf, rem)


class DiseaseSIR(DiseaseBase):
    def __init__(self, p, cmatrix, rng, fname, mode):
        super(DiseaseSIR, self).__init__(p, None, cmatrix, fname, mode)
        self.exposure = ExposureSEIRadd(p)
        rem = Removed(3)
        inf = Infected(2, DurationGenerator(
            p['infective_mean'] / (364.0 / p['t_per_year']), p['k'], rng))
        sus = Susceptible(0)
        self.add_states(sus, inf, rem)


class DiseaseSIRS(DiseaseBase):
    def __init__(self, p, cmatrix, rng, fname, mode):
        super(DiseaseSIRS, self).__init__(p, None, cmatrix, fname, mode)
        self.exposure = ExposureSEIRadd(p)
        rem = RemovedTemp(3, DurationGenerator(
            p['removed_mean'] / (364.0 / p['t_per_year']), p['k'], rng))
        inf = Infected(2, DurationGenerator(
            p['infective_mean'] / (364.0 / p['t_per_year']), p['k'], rng))
        sus = Susceptible(0)
        self.add_states(sus, inf, rem)

class DiseaseSIS(DiseaseBase):
    def __init__(self, p, cmatrix, rng, fname, mode):
        super(DiseaseSIS, self).__init__(p, None, cmatrix, fname, mode)
        self.exposure = ExposureSEIRadd(p)
        
        inf = Infected(2, DurationGenerator(
            p['infective_mean'] / (364.0 / p['t_per_year']), p['k'], rng))
        sus = Susceptible(0)
        self.add_states(sus, inf)

class DiseaseSIS_NoRecCases(DiseaseBase):
    def __init__(self, p, cmatrix, rng, fname, mode):
        super(DiseaseSIS_NoRecCases, self).__init__(p, None, cmatrix, fname, mode)
        self.exposure = ExposureSEIRadd(p)
        
        inf = Infected(2, DurationGeneratorNoRecCases(
            p['infective_mean'] / (364.0 / p['t_per_year']), p['k'], p['norecovery_cases_frac'], rng))
        sus = Susceptible(0)
        self.add_states(sus, inf)
    
class DiseaseMSIR(DiseaseBase):
    """
    A variant of standard SIR in which newborns enter a maternal protection class
    depending on the exposure history of their mother.
    
    maternal_inf_dur: maximum duration since natural infection that results in
        transmission of maternal antibodies
        
    maternal_vac_dur: maximum duration since vaccination that results in
        transmission of maternal antibodies
        
    We hypothesize that the latter is shorter than the former.
        
    As a greater proportion of mothers start to obtain immunity via vaccination
    than from natural infection, we would expect lower rates of maternal
    immunity and hence the potential for more disease among unvaccinated infants.
    """

    def __init__(self, p, cmatrix, rng, mat_inf_dur, mat_vac_dur, fname, mode):
        super(DiseaseMSIR, self).__init__(p, None, cmatrix, fname, mode)
        self.exposure = ExposureSEIRadd(p)
        rem = Removed(3)
        inf = Infected(2, DurationGenerator(
            p['infective_mean'] / (364.0 / p['t_per_year']), p['k'], rng))
        sus = Susceptible(1)
        mat = Maternal(0, DurationGenerator(
            p['maternal_mean'] / (364.0 / p['t_per_year']), p['k'], rng))
        self.add_states(mat, sus, inf, rem)

        self.mat_inf_dur = mat_inf_dur
        self.mat_vac_dur = mat_vac_dur

    def bd_update(self, t, births, deaths, imms):
        super(DiseaseMSIR, self).bd_update(t, births, deaths, imms)
        # test each new birth to see if they enter M rather than S

        # NB: current implementation (because it calls super function first) will
        # record in individual log that they entered S and then M on same timestep.
        # If this becomes an issue, may need to rewrite more of super function...
        for ind in births:
            mother = ind.parents[0]
            if (len(mother.infections) > 0 and (t - mother.infections[-1] < self.mat_inf_dur)
                or len(mother.vaccine_times) > 0 and (t - mother.vaccine_times[-1] < self.mat_vac_dur)):
                ind.next_state = self.states['M']
                self.tick(t, ind)


#class Susceptible_Maternal(Susceptible):
#    def __init__(self, successors, order, duration, maternal_dur):
#        super(Susceptible_Maternal, self).__init__(successors, order, duration)
#        self.maternal_dur = maternal_dur
#        
#    def enter(self, t, ind):
#        super(Susceptible_Maternal, self).enter(t, ind)
#        if ind.age==0 and ind.age_days==0:
#            mother = ind.parents[0]
#            # transfer a newborn into the maternal class if their mother's most recent
#            # exposure to infection was less than 'maternal_dur' timesteps ago
#            if len(mother.infections) > 0 and (t - mother.infections[-1] < self.maternal_dur):
#                ind.next_state = self.successors['M']


class DiseaseNetworkSEIR(DiseaseSEIR):
    def __init__(self, p, cmatrix, rng, fname, mode):
        super(DiseaseNetworkSEIR, self).__init__(p, cmatrix, rng, fname, mode)
        self.cnetwork = ContactNetwork()
        self.exposure = ExposureNetworkSEIRadd(p)
