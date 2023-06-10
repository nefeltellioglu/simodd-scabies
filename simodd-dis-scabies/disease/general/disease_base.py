"""
Template for updating disease state of a population.
"""
import tables as tb
from itertools import chain

import numpy as np


class DiseaseBase(object):
    def __init__(self, p, exposure, cmatrix, fname, mode,
                 basic_susceptible='S', basic_infection='I',
                 disease_states=tuple('I')):

        # logging: if True, store information on disease state transitions with each individual
        self.logging = p['logging']

        # halt: if True, terminate simulation when no individuals are in a disease state
        self.halt = p['halt']

        # external exposure rate: per person, per time step rate of disease introduction
        self.external_exposure_rate = p['external_exposure_rate']

        # exposure: function for calculating force of infection acting on an individual
        self.exposure = exposure

        # cmatrix / cnetwork: contact matrix (or network)
        self.cmatrix = cmatrix
        self.cnetwork = None

        # basic susceptible / infection; default states for (e.g.,) birth / importation
        self.basic_susceptible = basic_susceptible
        self.basic_infection = basic_infection

        # birth / death count: number of births / deaths in a timestep (WHY?)
        self.birth_count = 0
        self.death_count = 0

        # vaccines: set of vaccines to be applied during simulation (currently active or not)
        self.vaccines = {}

        # states: all disease model states
        self.states = {}

        # disease states: set of states that denote presence of disease (e.g., latent, infectious)
        self.disease_states = disease_states

        # infectious states: labels of states that contribute to force of infection
        self.infectious_states = None

        # by age: dictionary, keyed by state label, of number of individuals in that state per age class
        # for calculating community force of infection in an efficient fashion
        self.by_age = {}

        # h5file: HDF file for storing output of simulation
        self.h5file = tb.open_file(fname, mode)
        if mode in 'w':
            self.store_params(p)

        # obs_on: if True, record output on simulation
        self.obs_on = True

        # observers: data collection objects
        self.observers = {}

    ### Output file handling # # # # # # # # # # # #

    def store_params(self, p):
        """
        If not already present, store simulation parameters in output file and mark file as incomplete.
        """
        if '/params' not in self.h5file:
            self.h5file.create_group('/', 'params', 'Parameters')
        for k, v in list(p.items()):
            # rhdf5 has problems reading boolean values, therefore write as 0/1
            if type(v) is bool:
                self.h5file.set_node_attr('/params', k, 1 if v else 0)
            else:
                self.h5file.set_node_attr('/params', k, v)
        self.h5file.set_node_attr('/params', 'complete', 0)
        self.h5file.flush()

    def is_complete(self):
        """
        Check if current output file has been marked as complete.
        """
        return self.h5file.get_node_attr('/params', 'complete')

    def done(self, complete=False):
        """
        Close output file, marking as complete if specified.
        """
        if complete:
            self.h5file.set_node_attr('/params', 'complete', 1)
        self.h5file.close()

    def close(self):
        self.h5file.close()

    ### Initialisation and clean-up # # # # # # # # # # # #

    def set_counts(self, P):
        """
        Initialise state counts and I_by_age on basis of state of 
        individuals in population P.
        """
        # reset counts to zero
        for state in self.states.values():
            state.count = 0
        for label in self.infectious_states:
            if self.cmatrix:
                self.by_age[label] = [0] * len(self.cmatrix.age_classes)
            else:
                self.by_age[label] = [0]

        # set counts on basis of P
        for ind in P.I.values():
            # make sure individual is using this disease's state objects
            if ind.prev_state:
                ind.prev_state = self.states[ind.prev_state.label]
            ind.state = self.states[ind.state.label]
            if ind.next_state:
                ind.next_state = self.states[ind.next_state.label]
            # update state count
            self.states[ind.state.label].count += 1
            # update I_by_age dict

            if ind.state.label in self.infectious_states:
                self.by_age[ind.state.label][self.cmatrix.age_map[ind.age]] += 1
                ind.state.current.add(ind.ID)

        # update age-sorted lists of individuals
        self.cmatrix.init_age_classes(P)

    def add_states(self, *states):
        """Add a new disease state."""
        for state in states:
            self.states[state.label] = state
        self.__add_infectious_states()

    def __add_infectious_states(self):
        """Add infectious states."""
        self.infectious_states = [label for label in self.disease_states if self.states[label].infectious]
        for label in self.infectious_states:
            if self.cmatrix:
                self.by_age[label] = [0] * len(self.cmatrix.age_classes)
            else:
                self.by_age[label] = [0]

    def add_observers(self, *observers):
        """Add a new observer to the observers list."""
        for observer in observers:
            self.observers[observer.label] = observer

    def add_vaccines(self, *vaccines):
        """ Add new vaccines. """
        for vaccine in vaccines:
            self.vaccines[vaccine.label] = vaccine

    def seed_infection(self, t, P, cases, rng, seed_inds=None):
        """Seed initial infection (set everyone else to susceptible)."""
        for ind in P.I.values():
            ind.next_state = self.states[self.basic_susceptible]
            self.tick(t, ind)
        if not seed_inds:
            seed_inds = rng.sample(list(P.I.values()), cases)
        for ind in seed_inds:
            ind.next_state = self.states[self.basic_infection]
            ind.source = -3
            I_in_age, I_out_age = self.tick(t, ind)
            if I_in_age is not None and I_in_age >= 0:
                self.by_age[ind.state.label][self.cmatrix.age_map[I_in_age]] += 1

    ### Updating  # # # # # # # # # # # # # # # # # # # # #

    def bd_update(self, t, births, deaths, imms, rng):
        """
        Update state counts for births and deaths and immigration.
        
        Currently immigrants are treated as arriving susceptible.
        
        TODO: update this to handle, e.g., transfer of maternal immunity.
        """
        for ind in chain(births, imms):
            ind.next_state = self.states[self.basic_susceptible]
            self.tick(t, ind)
        for ind in deaths:
            if ind.state:
                ind.state.exit(t, ind)
                if ind.state.label in self.infectious_states:
                    self.by_age[ind.state.label][self.cmatrix.age_map[ind.age_at_infection]] -= 1
#            ind.state = None
        self.birth_count = len(births)
        self.death_count = len(deaths)

    def update(self, t, P, rng):
        """
        Update the disease state of the population.
        """
        self.check_vaccines(t, P, rng)
        # print "(update)", introduction
        
        cases = self.check_exposure(t, P, rng)
            
        introductions = self.external_exposure(t, P, rng)
        
        # print "(update)", introduction
        new_I = self.update_ind_states(t, P)

        self.update_observers(t, disease=self, pop=P,
                              cases=cases['infection'],
                              boosting=cases['boosting'],
                              introductions=introductions,
                              new_I=new_I, rng=rng)
        
            
            
        # if halting, halt if no exposed or infected individuals in population (eg for SIR)
        return self.halt and not sum([self.states[x].count for x in self.disease_states])

    def check_vaccines(self, t, P, rng):
        """
        Reset vaccine counts and check for successful vaccinations.
        """
        for v in self.vaccines.values():
            v.count = 0
            #v.vaccinate_pop(t, P, self.states, rng)
            if v.check_vaccination_period(t):
                v.vaccinate_pop(t, P, self.states, rng) #vaccination is applied at this time step
                if v.immediately_effective: #vaccine is immediately effective
                    #May42021 -added below line to update state of individuals-without updating the state counters-.
                    self.update_ind_states_after_vaccination(t, P)



    def check_exposure(self, t, P, rng):
        """
        Check all at risk individuals for potential exposure to infection.
        """

        # we will return a list of exposures (broken into infectious cases and boosting)
        cases = dict(infection=[], boosting=[])

        # return empty lists if there is currently no force of infection acting
        if sum([self.states[label].count for label in self.infectious_states]) == 0:
            return cases

        # create a set of individuals currently at risk
        # split depending upon whether network or matrix is being used to calculate community exposure
        pop_at_risk = set()
        comm = {}
        if self.cnetwork:
            # population at risk consists of household members and network neighbours of
            # people currently in an infectious state
            for cur_inf_state in self.infectious_states:
                for cur_I in self.states[cur_inf_state].current:
                    pop_at_risk.update([x for x in P.housemates(P.I[cur_I]) if x.state.at_risk])
                    pop_at_risk.update([x for x in self.cnetwork.get_contacts(P, P.I[cur_I]) if x.state.at_risk])
        else:
            # population at risk consists of potentially everybody
            pop_at_risk = [x for x in P.I.values() if x.state.at_risk]
            # compute exposure from community:
            # comm is the force of infection arising from infection in the community
            # EC[i] is a vector containing the contact rates between an individual in age group i
            # and individuals in each age group j, weighted by the (approximate) number of people
            # in age group j (as community mixing is density dependent).
            # I_by_age is the number of infected individuals in each age group j
            # total force of infection for each age class i is equal to the products of
            # weighted contact rate and the number of infected individuals, summed over
            # each age class i.
            for label in self.infectious_states:
                if self.cmatrix:
                    comm[label] = np.array([np.sum(self.cmatrix.EC[i] * self.by_age[label]) for i in range(101)])
                else:
                    # dummy line for non-age based mixing
                    comm[label] = np.array([float(self.states[label].count) / len(P.I) for _ in range(101)])

        # test exposure for each individual at risk
        for ind in pop_at_risk:
            # split depending upon whether network or matrix is being used to calculate community exposure
            if self.cnetwork:
                foi = self.exposure.calc_foi_fast(t, ind, P, self.cnetwork, rng)
            else:
                foi = self.exposure.calc_foi_fast(t, ind, P, comm, rng)
            exposure_type = ind.state.test_exposure(self.states, ind, foi, rng)
            #                    p, s = P.hh_parents_siblings(ind)
            #                    p_I = len([x for x in p if x.state.infectious])
            #                    s_I = len([x for x in s if x.state.infectious])
            if exposure_type == 'infection':
                ind.infections.append(t)
                cases['infection'].append(ind)
            elif exposure_type == 'boosting':
                cases['boosting'].append(ind)
        return cases

    def old_external_exposure(self, t, P, rng):
        """
        Check for external import of infection.
        """
        external_exposure_prob = self.external_exposure_rate * len(P.I)
        if external_exposure_prob < 0 and not sum([self.states[x].count for x in self.disease_states])\
                or rng.random() <= external_exposure_prob:
            tries = 0
            max_tries = 100
            while tries < max_tries:
                ind = rng.sample(list(P.I.values()), 1)[0]
                if ind.state.at_risk:
                    ind.next_state = self.states[self.basic_infection]
                    ind.source = -2
                    ind.infections.append(t)
                    return [ind]
                tries += 1
        return []
    
    def external_exposure(self, t, P, rng):
        """
        Check for external import of infection.
        """
        cases = []
        external_exposure_prob = self.external_exposure_rate * len(P.I)
        if external_exposure_prob < 0 and not sum([self.states[x].count for x in self.disease_states])\
                or rng.random() <= external_exposure_prob:
            tries = 0
            max_tries = 100
            while tries < max_tries:
                ind = rng.sample(list(P.I.values()), 1)[0]
                if ind.state.at_risk:
                    ind.next_state = self.states[self.basic_infection]
                    ind.source = -2
                    ind.infections.append(t)
                    cases.append(ind)
                    if len(cases) == 5:
                        return cases
                tries += 1
        return cases

    def update_ind_states(self, t, P):
        """
        Update disease status of all individuals.
        Returns a list of newly infectious ('symptomatic' individuals)
        """
        new_I = []
        # second loop updates current state
        for ind in P.I.values():
            old_state = ind.state.label
            I_in_age, I_out_age = self.tick(t, ind)
            new_state = ind.state.label

            if old_state not in self.infectious_states and new_state in self.infectious_states:
                new_I.append(ind)
            if not self.cmatrix: continue

            if I_in_age is not None and I_in_age >= 0 and new_state in self.infectious_states:
                self.by_age[new_state][self.cmatrix.age_map[I_in_age]] += 1
            if I_out_age is not None and I_out_age >= 0 and old_state in self.infectious_states:
                self.by_age[old_state][self.cmatrix.age_map[I_out_age]] -= 1

        return new_I
    #May42021 Nefel introduced update_ind_states_after vaccination function.
    def update_ind_states_after_vaccination(self, t, P):
        """
        Update disease status of all individuals after vaccination.
        Differences from update_ind_states():
        1.Does not return something. 
        2.Calls self.tick_after_vaccination(t, ind) instead of self.tick(t, ind)
        
        """
        new_I = []
        # second loop updates current state
        for ind in P.I.values():
            old_state = ind.state.label
            I_in_age, I_out_age = self.tick_after_vaccination(t, ind)
            new_state = ind.state.label

            if old_state not in self.infectious_states and new_state in self.infectious_states:
                new_I.append(ind)
            if not self.cmatrix: continue

            if I_in_age is not None and I_in_age >= 0 and new_state in self.infectious_states:
                self.by_age[new_state][self.cmatrix.age_map[I_in_age]] += 1
            if I_out_age is not None and I_out_age >= 0 and old_state in self.infectious_states:
                self.by_age[old_state][self.cmatrix.age_map[I_out_age]] -= 1

        #return new_I
    #May42021 Nefel introduced tick_after_vaccination function.
    def tick_after_vaccination(self, t, ind):
        """
        Update the state a single individual after vaccination.
        Differences from self.tick():
        1.Instead of calling ind.state.update(), it runs a 4 lines code chunk -where only keys which have negative
          values in the ind.counters are removed-.
        """
        recovered_age = None
        infected_age = None
        if ind.state != ind.next_state:
            if ind.state and ind.state.infectious:
                recovered_age = ind.age_at_infection
                ind.age_at_infection = None
            if ind.next_state and ind.next_state.infectious:
                infected_age = ind.age
                ind.age_at_infection = ind.age
            if self.logging:
                ind.add_log(t, ind.next_state.label,
                            'Individual entering state: %s' % ind.next_state.label)
            if ind.state:
                ind.state.exit(t, ind)
            ind.next_state.enter(t, ind)
        #difference from self.tick(): instead of calling ind.state.update() I run below 4 lines
        #ind.state.update(t, ind, self.states)
        labels = [key for key in ind.counters]
        for key in labels:
            if ind.counters[key] < 1:
                del ind.counters[key]
        return infected_age, recovered_age
    def tick(self, t, ind):
        """
        Update the state a single individual.
        
        Used to be a function of class individual, but wanted to be able to pass in self.states
        (so that state.successors could contain labels rather than references.)
        """
        recovered_age = None
        infected_age = None
        if ind.state != ind.next_state:
            if ind.state and ind.state.infectious:
                recovered_age = ind.age_at_infection
                ind.age_at_infection = None
            if ind.next_state and ind.next_state.infectious:
                infected_age = ind.age
                ind.age_at_infection = ind.age
            if self.logging:
                ind.add_log(t, ind.next_state.label,
                            'Individual entering state: %s' % ind.next_state.label)
            if ind.state:
                ind.state.exit(t, ind)
            ind.next_state.enter(t, ind)
        ind.state.update(t, ind, self.states)
        return infected_age, recovered_age

    def update_observers(self, t, **kwargs):
        """
        Store observed data (if observers are switched on).
        """
        if self.obs_on:
            for observer in self.observers.values():
                observer.update(t, **kwargs)

    ### Information  # # # # # # # # # # # # # # # # # # # # #

    def state_labels(self):
        """
        Return a list of state labels, in specified order.
        """
        return [v.label for v in sorted(
            self.states.values(), key=lambda x: x.order)]

    def state_colors(self):
        """
        Return a list of state labels, in specified order.
        """
        return [v.color for v in sorted(
            self.states.values(), key=lambda x: x.order)]

    def state_counts(self):
        """
        Return a list of state counts, in specified order.
        """
        return [v.count for v in sorted(
            self.states.values(), key=lambda x: x.order)]