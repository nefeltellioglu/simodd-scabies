"""
.. module:: simulation
.. moduleauthor:: Nic Geard <ngeard@unimelb.edu.au>
"""

import os
import tables as tb
from random import Random
from math import exp
from collections import defaultdict

from population.pop_hh import PopHH
from population.individual import Individual
from population.pop_gen import gen_hh_age_structured_pop, allocate_couples
from population.utils import sample_table, load_probs, \
    load_probs_new, load_age_rates, load_prob_tables, load_prob_list


def _adjust_prob(rate, t_per_year):
    """
    convert from an annual rate to a per-time-period probabiliy
    where time-period = 1/t_per_year

    :param rate: annual rate to convert.
    :type rate: double
    :param t_per_year: number of time periods per year.
    :type t_per_year: int
    """
    tp = 1.0 / t_per_year
    return 1 - pow(1 - rate, tp)


class Simulation(object):
    """
    Basic demographic simulation object.

    Handles updating of births, deaths, aging, immigration, and household structure (couple formation and
    disolution, leaving home, etc.)

    :param p: dictionary of simulation parameters.
    :type p: dict
    :param ind_type: the :class:`.individual.Individual` (sub)class stored by this population.
    :type ind_type: class
    :param create_pop: If `True` (default), create a random population; otherwise, this will need to be done later.
    :type create_pop: bool

    """

    def __init__(self, p, ind_type=Individual, create_pop=True, h5file=None):
        self.p = p  # local copy of parameters
        self.p_adj = {}  #
        self.rng = Random(self.p['seed'])  # random number generator
        self.ind_type = ind_type

        self.obs_on = True
        self.observers = {}

        self.growth_residue = 0.0

        # store this in a data/rate dictionary?
        self.hh_comp = None
        self.death_rates = None
        self.fertility_age_probs = None
        self.fertility_parity_probs = None
        self.dyn_years = None

        self._setup_params()

        self.h5file = h5file
        if self.h5file is None:
            self.h5file = tb.open_file(os.path.join(p['prefix'], 'population.hd5'), 'w')
            self._store_params(p)

        self.P = None
        if create_pop:
            self.create_population()

        # map of times to mothers due
        self.preg_schedule = defaultdict(list)

    def _setup_params(self):
        """
        Reset various simulation parameters.

        (could be cleaner what is going on here between initialization
        and resetting...)
        """
        self.growth_residue = 0.0
        self._load_demographic_data()

    def _store_params(self, p):
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

    def create_population(self):
        """
        Create a population according to specified age and household size
        distributions.
        """
        self._setup_params()
        self.P = PopHH(self.ind_type, self.p['logging'])
        gen_hh_age_structured_pop(self.P, self.p['pop_size'], self.hh_comp,
                                  self.age_dist, self.p['age_cutoffs'], self.rng)
        allocate_couples(self.P)

    def add_observers(self, *observers):
        """Add a new observer to the observers list."""
        for observer in observers:
            self.observers[observer.label] = observer

    @staticmethod
    def _parse_age_rates(filename, factor, final):
        """
        Parse an age-year-rate table to produce a dictionary, keyed by age, 
        with each entry being a list of annual rates (by year).
        
        Setting final to 'True' appends an age 100 rate of >1 (e.g., to 
        ensure everyone dies!
        """

        dat = load_age_rates(filename)
        rates = {}
        for line in dat:
            rates[line[0]] = [x * factor for x in line[1:]]
        if final:
            rates[100] = [100 for _ in dat[0][1:]]  # added to try and fix appearance of individuals >100y
            rates[101] = [100 for _ in dat[0][1:]]  # everybody dies...
        return rates

    def _load_demographic_data(self):
        """
        Load data on age-specific demographic processes (mortality/fertility)
        and adjust event probabilities according to time-step.

        All paths, parameter values, etc. are contained in the parameter dictionary
        passed in when the :class:`Simulation` object is created.
        """

        # load household size distribution and age distribution
        self.hh_comp = load_probs(os.path.join(self.p['resource_prefix'],
                                               self.p['hh_composition']), False)
        self.p['age_cutoffs'] = [int(x) for x in self.hh_comp[0][1:][0]]  # yuk!
        self.age_dist = load_probs(os.path.join(self.p['resource_prefix'],
                                                self.p['age_distribution']))

        annual_factor = 1.0 / self.p['t_per_year']

        # load and scale MORTALITY rates
        self.death_rates = {
            0: self._parse_age_rates(os.path.join(
                self.p['resource_prefix'],
                self.p['death_rates_m']), annual_factor, True),
            1: self._parse_age_rates(os.path.join(
                self.p['resource_prefix'],
                self.p['death_rates_f']), annual_factor, True)}

        # load FERTILITY age probs (don't require scaling) for closed pops
        self.fertility_age_probs = load_prob_tables(os.path.join(
            self.p['resource_prefix'],
            self.p['fertility_age_probs']))
        self.fertility_parity_probs = load_probs_new(os.path.join(
            self.p['resource_prefix'],
            self.p['fertility_parity_probs']))

        # load and scale leave/couple/divorce and growth rates
        if self.p['dyn_rates']:
            # rates will be a list of annual values
            self.p['leaving_probs'] = load_prob_list(os.path.join(
                self.p['resource_prefix'], self.p['leaving_prob_file']))
            self.p['couple_probs'] = load_prob_list(os.path.join(
                self.p['resource_prefix'], self.p['couple_prob_file']))
            self.p['divorce_probs'] = load_prob_list(os.path.join(
                self.p['resource_prefix'], self.p['divorce_prob_file']))
            self.p['growth_rates'] = load_prob_list(os.path.join(
                self.p['resource_prefix'], self.p['growth_rate_file']))
            self.p['imm_rates'] = load_prob_list(os.path.join(
                self.p['resource_prefix'], self.p['imm_rate_file']))

            # print self.p['growth_rates'][0]

            self.p_adj['leaving_probs'] = [_adjust_prob(x, self.p['t_per_year'])
                                           for x in self.p['leaving_probs']]
            self.p_adj['couple_probs'] = [_adjust_prob(x, self.p['t_per_year'])
                                          for x in self.p['couple_probs']]
            self.p_adj['divorce_probs'] = [_adjust_prob(x, self.p['t_per_year'])
                                           for x in self.p['divorce_probs']]
            self.p_adj['growth_rates'] = [_adjust_prob(x, self.p['t_per_year'])
                                          for x in self.p['growth_rates']]
            self.p_adj['imm_rates'] = [_adjust_prob(x, self.p['t_per_year'])
                                       for x in self.p['imm_rates']]

            # print self.p['growth_rates'][0]

            self.dyn_years = min(len(self.death_rates[0][0]) - 1,
                                 len(self.fertility_age_probs) - 1,
                                 len(self.p_adj['leaving_probs']) - 1,
                                 len(self.p_adj['couple_probs']) - 1,
                                 len(self.p_adj['divorce_probs']) - 1,
                                 len(self.p_adj['growth_rates']) - 1)

        else:
            # adjust demographic event probabilities according to time step
            self.p_adj['couple_probs'] = [_adjust_prob(
                self.p['couple_prob'], self.p['t_per_year'])]
            self.p_adj['leaving_probs'] = [_adjust_prob(
                self.p['leaving_prob'], self.p['t_per_year'])]
            self.p_adj['divorce_probs'] = [_adjust_prob(
                self.p['divorce_prob'], self.p['t_per_year'])]
            self.p_adj['growth_rates'] = [_adjust_prob(
                self.p['growth_rate'], self.p['t_per_year'])]
            self.p_adj['imm_rates'] = [_adjust_prob(
                self.p['imm_rate'], self.p['t_per_year'])]

    def _update_individual_demo(self, t, ind, index=0):
        """
        Update individual ind; check for death, couple formation, leaving home
        or divorce, as possible and appropriate.
        """

        death = None
        birth = None

        # DEATH / BIRTH:
        if self.rng.random() > exp(-self.death_rates[ind.sex][ind.age][index]):
            # currently pregnant women are 'immune' from dying
            if ind in self.P.preg_current:
                return death, birth
            death = ind
            mother = ind
            # make sure the dead or currently pregnant individual isn't selected as mother
            while mother is ind or mother in self.P.preg_current:
                mother = self._choose_mother(index)
            if mother == "error":
                return "error", "error"

            # trigger death and reallocate any orphan children
            orphans = self.P.death(t, ind)
            for cur_dep in orphans:
                if cur_dep.age > self.p['leaving_age']:
                    self.P.leave_home(t, cur_dep)
                else:
                    hh = self._choose_household(cur_dep)
                    self.P.allocate_orphan(t, cur_dep, hh)

            if self.p['preg']:
                # add new mother to pregnancy schedule
                due_time = t + int(280.0 / (364.0 / self.p['t_per_year']))
                self.P.preg_current[mother] = due_time
                self.preg_schedule[due_time].append(mother)
            else:
                # birth occurs immediately
                birth = self.P.birth(t, mother, 0 if self.rng.random() < 0.5 else 1)

        # COUPLE FORMATION:
        elif self.p['couple_age'] < ind.age < self.p['couple_age_max'] \
                and not ind.partner \
                and self.rng.random() < self.p_adj['couple_probs'][index]:
            partner = self._choose_partner(ind)
            if partner:
                self.P.form_couple(t, ind, partner)

        # LEAVING HOME:
        elif ind.age > self.p['leaving_age'] \
                and ind.with_parents \
                and not ind.partner \
                and self.rng.random() < self.p_adj['leaving_probs'][index]:
            self.P.leave_home(t, ind)

        # DIVORCE:
        elif self.p['divorce_age'] < ind.age < self.p['divorce_age_max'] \
                and ind.partner \
                and self.rng.random() < self.p_adj['divorce_probs'][index]:
            self.P.separate_couple(t, ind)

        # ELSE: individual has a quiet year...
        return death, birth

    def _choose_mother(self, index):
        """
        Choose a new mother on the basis of fertility rates.

        :param index: the index of the current rate set to use (for dynamic rates)
        :type index: int

        ..note::

            There is still a very small possibility (in *very* small
            populations) that this will fail due to candidates remaining
            forever empty.  This is currently handled by propagating an
            error code back up, that can be used to trigger a re-attempt.

        """

        candidates = []
        attempts = 0
        max_attempts = 500  # before restarting with a new population
        while not candidates:
            tgt_age = int(sample_table(
                self.fertility_age_probs[index], self.rng)[0])
            tgt_prev_min = 0
            tgt_prev_max = 100
            if self.p['use_parity']:
                tgt_prev_min = int(sample_table(
                    self.fertility_parity_probs[(tgt_age-15)//5], self.rng)[0])
                # effectively transform 5 into 5+
                tgt_prev_max = tgt_prev_min if tgt_prev_min < 5 else 20
            tgt_set = self.P.individuals_by_age(tgt_age, tgt_age)
            # print len(tgt_set)
            # print [(x.sex, x.can_birth(), not x.with_parents, len(x.children)) for x in tgt_set]
            candidates = [
                x for x in tgt_set
                if x.sex == 1 and x.can_birth() and not x.with_parents
                and tgt_prev_min <= len(x.children) <= tgt_prev_max
                and x not in self.P.preg_current
            ]
            # print candidates
            attempts += 1
            if attempts >= max_attempts:
                # probably a better way to do this, but this "error" is
                # currently propagated upwards (eventually to sim_epi)
                return "error"
        return self.rng.choice(candidates)

    def _choose_partner(self, ind):
        """
        Choose a partner for i_id, subject to parameter constraints.

        :param ind: the first partner in the couple.
        :type ind: ind_type
        :returns: partner if successful, otherwise None.
        """

        mean_age = ind.age + self.p['partner_age_diff'] \
            if ind.sex == 0 else ind.age - self.p['partner_age_diff']
        tgt_age = 0
        candidates = []
        while tgt_age < self.p['min_partner_age']:
            tgt_age = int(self.rng.gauss(mean_age, self.p['partner_age_sd']))
            tgt_set = self.P.individuals_by_age(tgt_age, tgt_age)
            candidates = [
                x for x in tgt_set
                if not x.partner and x.sex != ind.sex
                and x not in self.P.groups['household'][ind.groups['household']]
            ]

        # abort if no eligible partner exists
        return None if not candidates else self.rng.choice(candidates)

    def _choose_household(self, ind):
        """
        Process orphans who result when the last remaining adult guardian
        in their household dies.  If they are above 'adult-age' cutoff, place
        them in a new single household, otherwise, reallocate them to an
        existing family household (with at least one other child).
        """
        cur_hh = ind.groups['household']
        # choose destination household from among family households
        candidates = self.P.groups_by_min_size('household', 3)
        tgt_hh = cur_hh
        while tgt_hh == cur_hh:
            tgt_hh = self.rng.sample(candidates, 1)[0]
        return tgt_hh

    def update_all_demo(self, t):
        """
        Carry out a single update of all demographic aspects population.

        :param t: the current time step.
        :type t: int
        :returns: a tuple containing lists of births, deaths, immigrants and birthdays

        """

        # age each individual by appropriate number of days
        birthdays = self.P.age_population(364 // self.p['t_per_year'])

        deaths = []
        births = []

        # calculate index for fertility and mortality rates
        # basically: use first entry for burn-in, then one entry every 
        # 'period' years, then use the final entry for any remaining years.
        index = min(max(
            0, (t - (self.p['demo_burn'] * self.p['t_per_year'])) //
               (self.p['t_per_year'])), self.dyn_years) \
            if self.p['dyn_rates'] else 0

        # print cur_t / self.p['t_per_year'], \
        #     index, self.p_adj['growth_rates'][index], \
        #     len(self.P.I)

        cur_inds = list(self.P.I.values())
        for ind in cur_inds:
            death, birth = self._update_individual_demo(t, ind, index)
            if death == "error" and birth == "error":
                return "error", "error", "error", "error"
            if death:
                deaths.append(death)
            if birth:
                births.append(birth)

        # trigger delayed (due to pregnancy births)
        for mother in self.preg_schedule[t]:
            # create new individuals
            new_ind = self.P.birth(t, mother, 0 if self.rng.random() < 0.5 else 1)
            births.append(new_ind)
            # remove mother from pregnancy list
            del self.P.preg_current[mother]
        del self.preg_schedule[t]

        # population growth
        ####################

        # print "growth rate (%d) = %f: %d" % (index, self.p_adj['growth_rates'][index], len(self.P.I))

        # calculate number of new individuals to add (whole and fraction)
        new_individuals = len(self.P.I) * self.p_adj['growth_rates'][index]
        # get whole part of new individuals
        new_now = int(new_individuals)
        # add fractional part to residue accumulation
        self.growth_residue += (new_individuals - new_now)
        # grab any new 'whole' individuals
        new_residue = int(self.growth_residue)
        new_now += new_residue
        self.growth_residue -= new_residue

        # create the new individuals
        for _ in range(new_now):
            mother = self._choose_mother(index)
            if mother == "error":
                return "error", "error", "error", "error"
            new_ind = self.P.birth(t, mother, 0 if self.rng.random() < 0.5 else 1)
            births.append(new_ind)

        # immigration
        ##############

        imm_count = 0
        imm_tgt = int(len(self.P.I) * self.p_adj['imm_rates'][index])
        source_hh_ids = []
        immigrants = []
        while imm_count < imm_tgt:
            hh_id = self.rng.choice(list(self.P.groups['household'].keys()))
            imm_count += len(self.P.groups['household'][hh_id])
            source_hh_ids.append(hh_id)
        for hh_id in source_hh_ids:
            new_hh_id = self.P.duplicate_household(t, hh_id)
            immigrants.extend(self.P.groups['household'][new_hh_id])

        return births, deaths, immigrants, birthdays

    def update_observers(self, t, **kwargs):
        """
        Store observed data (if observers are switched on).
        """
        if self.obs_on:
            for observer in self.observers.values():
                observer.update(t, **kwargs)

    def done(self, complete=False):
        """
        Close output file, marking as complete if specified.
        """
        if complete:
            self.h5file.set_node_attr('/params', 'complete', 1)
        self.h5file.close()
