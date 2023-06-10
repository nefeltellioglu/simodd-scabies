"""
Observer for information on every agent.
"""
from itertools import chain
import os
from math import ceil
import tables as tb
from collections import defaultdict
from scipy.optimize import curve_fit

from pylab import NaN

from population.pop_info import hh_size, age_dist, age_dist_month
from disease.observers.obs_base import Observer
from disease.experiments.output_disease import *

month_days = [0, 30, 61, 89, 120, 150, 181, 211, 242, 273, 303, 334]


class DenomCase(tb.IsDescription):
    """
    A 'denominator' case record.
    """

    age = tb.UInt8Col()
    age_days = tb.UInt16Col()
    maternal_vacc = tb.UInt8Col()
    cocoon_vacc = tb.UInt8Col()
    birth_order = tb.UInt8Col()


class Case(tb.IsDescription): #SIS
    """
    An individual record. It doesn't have to be an infected case.
    """
#zzz
    #t = tb.UInt32Col()
    #t_E = tb.UInt32Col()
    #t_R = tb.UInt32Col()
    #dur_E = tb.UInt32Col()
    #dur_I = tb.UInt32Col()
    ID = tb.UInt32Col()
    state = tb.StringCol(4)
    age = tb.UInt8Col()
    sex = tb.UInt8Col()
    #age_days = tb.UInt16Col()
    hh_size = tb.UInt8Col()
    #hh_age = tb.UInt8Col()
    #hh_type = tb.StringCol(16)
    #num_siblings = tb.UInt8Col()
    #num_parents = tb.UInt8Col()
    birth_order = tb.UInt8Col()
    first = tb.Int8Col()
    #hh_frac = tb.Float32Col()
    hh_id = tb.UInt32Col()
    hh_source = tb.Int8Col()
    #ext_contacts = tb.UInt8Col()
    #ext_infections = tb.UInt8Col()
    #hh_infections = tb.UInt8Col()
    source_id = tb.Int32Col()
    source_age = tb.UInt8Col()
    #hh_at_risk = tb.UInt32Col()  # number of susceptible individuals upon hh introduction (ie, hh_frac == 0)
    hh_index_case = tb.Int8Col()  # True if there are no infectious individuals in household at time of case
    hh_index_case_first = tb.Int8Col()  # True if there have *never* been any infectious individuals in household
    maternal_vacc = tb.UInt8Col()  # True if mother received maternal vaccination
    #cocoon_vacc = tb.UInt8Col()  # True if mother received cocoon vaccination
    primary_vacc = tb.UInt8Col()  # True if ind receiving primary vaccination course


class AllObserver(Observer): #SIS
    """
    An observer for recording all occurring cases.
    """
    def __init__(self, h5file, t_per_year=364, 
                 #entry_state='I', #SIS
                 store_recovery=False, store_source=True):
        self.age_dists = None
        self.age_dist = None
        self.age_dists_month = None
        self.age_dist_month = None
        self.age_dists_month_first = None
        self.age_dists_month_subs = None
        self.age_dist_month_first = None
        self.age_dist_month_subs = None
        self.hh_size_dists = None
        self.hh_size_dist = None
        self.pop_dists_group = None
        self.t_per_year = t_per_year
        #self.entry_state = entry_state
        self.store_recovery = store_recovery
        self.store_source = store_source
        if self.store_recovery:
            self.open_cases = defaultdict(list)  # open cases
        super(AllObserver, self).__init__(h5file, 'cases', Case, 'Case Data')

    def create_storage(self, description, title):
        """
        Called when initialised with a h5file opened in 'w'rite mode.
        Adds tables for storage of data (plus local links)
        """

        Observer.create_storage(self, description, title)
        group = self.h5file.get_node('/', 'cases')
        self.age_dists = self.h5file.create_table(group, 'age_dists',
                                                 {'t': tb.UInt32Col(), 'dist': tb.UInt32Col(shape=(101))},
                                                 'Age Distributions')
        self.age_dist = self.age_dists.row
        self.age_dists_month = self.h5file.create_table(
            group, 'age_dists_month',
            {'t': tb.UInt32Col(), 'dist': tb.UInt32Col(shape=(60))}, 'Age Distributions (months)')
        self.age_dist_month = self.age_dists_month.row
        self.age_dists_month_first = self.h5file.create_table(
            group, 'age_dists_month_first',
            {'t': tb.UInt32Col(), 'dist': tb.UInt32Col(shape=(12))}, 'Age Distributions (months, first child)')
        self.age_dist_month_first = self.age_dists_month_first.row
        self.age_dists_month_subs = self.h5file.create_table(
            group, 'age_dists_month_subs',
            {'t': tb.UInt32Col(), 'dist': tb.UInt32Col(shape=(12))}, 'Age Distributions (months, subsequent child)')
        self.age_dist_month_subs = self.age_dists_month_subs.row
        self.hh_size_dists = self.h5file.create_table(
            group, 'hh_size_dists',
            {'t': tb.UInt32Col(), 'dist': tb.UInt32Col(shape=(20))}, 'Household Size Distributions')
        self.hh_size_dist = self.hh_size_dists.row

        self.pop_dists_group = self.h5file.create_group(group, 'pop_dists', 'Denominator populations')


    def load_storage(self):
        """
        Called when initialised with a h5file opened in 'r'ead or 'a'ppend mode.
        Creates local links to existing tables.
        """

        Observer.load_storage(self)
        self.age_dist = self.h5file.root.cases.age_dists.row
        self.hh_size_dist = self.h5file.root.cases.hh_size_dists.row
        self.age_dists = self.h5file.root.cases.age_dists
        self.age_dist = self.age_dists.row
        #if self.h5file.root.cases.__contains__('age_dists_month'):
        self.age_dists_month = self.h5file.root.cases.age_dists_month
        self.age_dist_month = self.age_dists_month.row
        self.hh_size_dists = self.h5file.root.cases.hh_size_dists
        self.hh_size_dist = self.hh_size_dists.row

    ### NB: following won't work, as create_storage needs arguments.
    ### If this is actually used anywhere, will need fixing!!!
    ### 10/1/14: looks ok to me? only used in estimate_R.py
    def _reset(self):
        """
        Remove all data and re-create empty storage.
        """
        self.h5file.remove_node('/', 'cases', recursive=True)
        self.create_storage(Case, 'Case Data')

    def _add_new_pop_dist(self, t, pop):
        filters = tb.Filters(complevel=9)
        cur_pop = self.h5file.create_table(self.pop_dists_group, 't%d'%t, DenomCase, filters=filters)
        row = cur_pop.row
        for ind in pop.individuals_by_age(0):
            row['age'] = ind.age
            row['age_days'] = ind.age_days
            row['birth_order'] = ind.birth_order
            row['maternal_vacc'] = 1 if (ind.parents and 'maternal' in ind.parents[0].vaccines_received) else 0
            row['cocoon_vacc'] = 1 if (ind.parents and 'cocoon' in ind.parents[0].vaccines_received) else 0
            row.append()
        cur_pop.flush()

    def _add_new_case(self, t, disease, pop, ind):
        """
        Append details of a new case to data table.
        """
        #zzz
        # p, s = pop.hh_parents_siblings(ind)

        #print 'OOO', ind.durations.get('E', 0), ind.durations.get('I', 0)

        #self.row['t_E'] = t
        #self.row['t'] = t + ind.durations.get('E', 0)
        if self.store_recovery:
            self.row['t_R'] = t + ind.durations.get('E', 0) + ind.durations.get('I', 0)
            self.row['dur_E'] = ind.durations.get('E', 0)
            self.row['dur_I'] = ind.durations.get('I', 0)
        self.row['ID'] = ind.ID
        self.row['state'] = ind.state.label
        self.row['age'] = ind.age
        #self.row['age_days'] = ind.age_days
        self.row['sex'] = ind.sex
        #self.row['age_days'] = ind.age_days
        self.row['hh_size'] = hh_size(pop, ind)
        #self.row['hh_age'] = pop.hh_age(t, ind)
        #self.row['hh_type'] = pop.get_hh_type(ind)
        #self.row['num_siblings'] = s
        #self.row['num_parents'] = p
        self.row['birth_order'] = ind.birth_order
        self.row['first'] = 1 if (len(ind.infections) == 1) else 0
        self.row['hh_id'] = ind.groups['household']
        #self.row['hh_frac'] = ind.hh_frac
        self.row['hh_source'] = 1 if ind.hh_source else 0
        #if disease.cnetwork:
        #    self.row['ext_contacts'] = disease.cnetwork.G.degree(ind.ID)
        #self.row['ext_infections'] = 0
        #self.row['hh_infections'] = 0
        self.row['source_id'] = ind.source
        self.row['source_age'] = pop.I[ind.source].age if ind.source > 0 else 0

        # if we know that individual was definitely infected from outside of the household, store the number
        # of other people in their household who were susceptible at that point in time. Add one for infected individual
        #self.row['hh_at_risk'] = len(
        #    [x.ID for x in pop.housemates(ind) if x.state.at_risk]) + 1 if ind.hh_frac <= 0.1 else 0

        housemates = pop.groups['household'][ind.groups['household']]

        self.row['hh_index_case'] = 1 if (not ind.hh_source and
                                          (sum(1 for x in housemates
                                               if x.state.infectious) == 0)) else 0

        self.row['hh_index_case_first'] = 1 if (not ind.hh_source and
                                                (sum(1 for x in housemates
                                                      if len(x.infections) > 0) == 0)) else 0

        self.row['maternal_vacc'] = 1 if (ind.parents and 'maternal' in ind.parents[0].vaccines_received) else 0
        #self.row['cocoon_vacc'] = 1 if (ind.parents and 'cocoon' in ind.parents[0].vaccines_received) else 0
        self.row['primary_vacc'] = 1 if ('primary_2m' in ind.vaccines_received) else 0

        self.row.append()
        self.data.flush()

        #if self.store_source:
        #    # update ext_infections count for source of infection
        #    s_id = ind.source
        #    if not s_id or s_id < 0:
        #        return  # no valid source
        #    if tb.__version__[0] == '3':
        #        source_row = self.data.get_where_list('ID == s_id')
        #    else:
        #        source_row = self.data.get_where_list('ID == s_id')
        #    if len(source_row) >= 1:
        #        if ind.hh_source:
        #            self.data.cols.hh_infections[source_row[-1]] += 1
        #        else:
        #            self.data.cols.ext_infections[source_row[-1]] += 1

    def update(self, t, disease, pop, cases, introductions=None, **kwargs):
        """
        Update details of all new cases occurring this time step.
        Also (optionally) store additional population data.

        :param t: The current time step
        :type t: int
        :param disease: The disease
        :type t: :class:`DiseaseBase`
        :param pop: The population
        :type pop: :class:`PopBase`
        :param cases: A list of individuals infected during this time step
        :type cases: list
        """

        # add cases (and introductions too, if any were passed)
        # print "(obs update)", introduction
        # print "# intro =", len([introduction] if introduction else [])
        #for ind in (cases + introductions) if introductions else cases: #SIS
        for ind in pop.I.values(): #SIS
            self._add_new_case(t, disease, pop, ind)
            if self.store_recovery:
                # build a list of individuals who will need their infectious durations updated when available
                self.open_cases[t + ind.durations.get('E', 0) + 1].append(ind.ID)
        # add age and household size distributions (4 times per year)
        if t % self.t_per_year in range(0, self.t_per_year, self.t_per_year // 4):
            # pop dists used for denominators when calculating incidence rates
#            self._add_new_pop_dist(t, pop)
            self.age_dist['t'] = t
            self.age_dist['dist'] = age_dist(pop.I, norm=False)[0]
            self.age_dist.append()
            self.age_dist_month['t'] = t
#            self.age_dist_month['dist'] = pop.age_dist_month(bin_days=91, max_age=5, norm=False)[0]
            self.age_dist_month['dist'] = age_dist_month(pop.I,
                                                         bin_days=month_days,
                                                         max_age=5, norm=False)[0]
            self.age_dist_month.append()
            self.hh_size_dist['t'] = t
            #            self.hh_size_dist['dist'] = pop.group_size_dist('household', 20, norm=False)[0]
            I_by_hh_size = pop.all_individuals_by_group_size('household', 20)
            self.hh_size_dist['dist'] = [len(I_by_hh_size[x]) for x in sorted(I_by_hh_size)]
            self.hh_size_dist.append()

        if self.store_recovery:
            # 're'-process individuals moving from 'E' to 'I' to add in their infectious duration
            # ('t_R') time entering R
            for cur_id in self.open_cases[t]:
                source_row = self.data.get_where_list('ID == %s' % cur_id)[-1]
                cur_ind = pop.I.get(cur_id, None)
                if not cur_ind:
                    cur_ind = pop.graveyard.get(cur_id, None)
                if not cur_ind:
                    print(("Error: individual %d not in I or graveyard..." % cur_id))
                    exit()
                self.data.cols.t_R[source_row] = self.data.cols.t[source_row] + cur_ind.durations.get('I', 0)
                self.data.cols.dur_I[source_row] = cur_ind.durations.get('I', 0)

            # remove current time step's open cases
            del self.open_cases[t]

        # Need to flush occasionally to avoid PerformanceWarnings
        self.h5file.flush()

    #########################################################################
    ### query and analysis functions ###

    def filter_data_by_state(self, st=None, et=None):
        """
        Returns a dictionary mapping state name (e.g., I, Ie, etc.) to data
        in that state.  Used if there is more than one type of infectious
        state.
        
        NB: note that values in filtered are iterators.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :returns: A dictionary mapping state names to case iterators

        """

        uniqvals = set()
        for x in self.data:
            uniqvals.add(x['state'])
        filtered = {}
        for cur_state in uniqvals:
            if st is None or et is None:
                filtered[cur_state] = self.data.where('state == cur_state')
            else:
                filtered[cur_state] = self.data.where('(state == cur_state) & (t>=st) & (t<et)')
        return filtered

    def get_num_cases(self, st=None, et=None):
        """
        Get number of cases in specified period (all by default).

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :returns: the number of cases in the specified time period

        """

        if st is None or et is None:
            return self.data.nrows
        else:
            return len([x for x in self.data.where('(t>=st) & (t<et)')])

    def get_case_IDs(self, st=None, et=None):
        """
        Get IDs of cases in specified period (all by default).

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :returns: a list of case IDs.
        """

        if st is None or et is None:
            return self.data.col('ID')
        else:
            return self.data.col('ID').where('(t>=st) & (t<et)')

    def get_new_data_by_time(self, dt=1, loc=None):
        """
        Returns a dictionary mapping times to the number of new infections occurring at that time.

        :param dt: time step increment
        :type dt: int
        :param loc: optionally, filter by location ('hh' or 'comm')
        :type: string
        :returns: A dictionary mapping timesteps to case counts.
        """
        case_t = defaultdict(int)
        if not loc:
            for x in self.data:
                case_t[(x['t'] // dt)] += 1
        elif loc == 'hh':
            for x in self.data.where('hh_source==1'):
                case_t[(x['t'] // dt)] += 1
        elif loc == 'comm':
            for x in self.data.where('hh_source==0'):
                case_t[(x['t'] // dt)] += 1
        return case_t

    def get_new_data_by_time_first(self, dt=1):
        """
        Returns a dictionary mapping times to the number of new *first* infections occurring at that time.

        :param dt: time step increment
        :type dt: int
        :returns: A dictionary mapping timesteps to case counts.
        """
        case_t = defaultdict(int)
        for x in self.data.where('first==1'):
            case_t[(x['t'] // dt)] += 1
        return case_t

    def get_new_case_counts_by_time(self, st, et, dt, loc=None):
        """
        Returns a list of case counts per (arbitrary) time unit, optionally in
        specific locations.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :param dt: time step increment
        :type dt: int
        :param loc: optionally, filter by location ('hh' or 'comm')
        :type: string
        :return: A list of case counts.
        """
        case_t = self.get_new_data_by_time(dt, loc)
        new_data = []
        upper = et // dt if et % dt == 0 else et // dt+1
        for cur_t in range(st // dt, upper):
            new_data.append(case_t.get(cur_t, 0))
        return new_data

    def get_new_case_counts_by_time_first(self, st, et, dt):
        """
        Returns a list of *first* case counts per (arbitrary) time unit, optionally in
        specific locations.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :param dt: time step increment
        :type dt: int
        :param loc: optionally, filter by location ('hh' or 'comm')
        :type: string
        :return: A list of case counts.
        """
        case_t = self.get_new_data_by_time_first(dt)
        new_data = []
        upper = et // dt if et % dt == 0 else et // dt+1
        for cur_t in range(st // dt, upper):
            new_data.append(case_t.get(cur_t, 0))
        return new_data

    def get_cum_case_count(self, st, et, dt=364):
        """
        Return a list of the cumulative number of cases occurring in each time period.

        :param et: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :param dt: time step increment
        :type dt: int
        :return: a list of cumulative case counts
        """

        return [self.get_num_cases(st, x + 1) for x in range(st, et, dt)]  # / 364))]

    def get_first_case_time(self):
        """
        Get time of first case.

        :returns: the time of the first case occurring.
        """

        return self.data[0]['t_E']

    def get_series(self, func, st=None, et=None, dt=1):
        """
        Get a time-series of the values generated by func.
        :param func: A function that takes a start and end time, and returns some value calculated over this period.
        :type func: function pointer
        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :param dt: time step increment
        :type dt: int
        :return: A list of tuples, containing mean/SD age of first infection for that time period.

        """

        cur_st = st
        cur_et = min(cur_st + dt, et)
        series = []
        while cur_st < et:
            series.append(func(cur_st, cur_et))
            cur_st = cur_et
            cur_et += dt
        return series

    def get_series_dict(self, func, st=None, et=None, dt=1):
        cur_st = st
        cur_et = min(cur_st + dt, et)
        series = defaultdict(list)
        while cur_st < et:
            cur_data = func(cur_st, cur_et)
            for k, v in list(cur_data.items()):
                series[k].append(v)
            cur_st = cur_et
            cur_et += dt
        return series


    #########################################################################
    ## age at first infection 

    def get_avg_age_first_infection(self, st=None, et=None):
        """
        Get mean/SD of ages of first infection in specified time period (all by default).

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :returns: a tuple containing (mean/SD) age of first infection
        """

        if st is None or et is None:
            ages = [x['age'] for x in self.data.where('first')]
        else:
            ages = [x['age'] for x in self.data.where('first & (t>=st) & (t<et)')]
        return np.mean(ages), np.std(ages) if len(ages) > 0 else (np.nan, np.nan)

    def get_hh_avg_age_first_infection(self, st=None, et=None):
        """
        Get a list of ages of first infection in specified time period,
        broken down by household size.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :returns: a list of ages of first infection (by hh size) for cases in the specified time period
        """

        # collection lists of ages
        ages = defaultdict(list)
        size_age = [(x['hh_size'], x['age']) for x in self.data.where('first')] \
            if (st is None or et is None) \
            else [(x['hh_size'], x['age']) for x in self.data.where('first & (t>=st) & (t<et)')]
        for cur_case in size_age:
            # store ages for hh of size (1..6, 7+)
            ages[min(cur_case[0], 7)].append(cur_case[1])
        # calculate average ages
        avg_ages = defaultdict(list)
        for cur_hh_size in range(1, 8):
            avg_ages[cur_hh_size].append(
                NaN if cur_hh_size not in ages
                else np.mean(ages[cur_hh_size]))
        return avg_ages

    #########################################################################
    ## fraction of infections attributable to households

    def get_avg_hh_frac(self, st=None, et=None):
        """
        Get avg fraction of infection attributable to hh transmission.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :returns: a tuple containing (mean/SD) fraction of infection attributable to hh transmission
        """

        if st is None or et is None:
            fracs = [x['hh_frac'] for x in self.data]
        else:
            fracs = [x['hh_frac'] for x in self.data.where('(t>=st) & (t<et)')]
        return (np.mean(fracs), np.std(fracs)) if len(fracs) > 0 else (np.nan, np.nan)

    def get_hh_hh_frac(self, st=None, et=None):
        """
        Get a household fractions in specified time period,
        broken down by household size.
        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :returns: a list of average fractions of infection attributable to hh transmission by hh size
        """
        hh_fracs = defaultdict(list)
        size_hh_frac = [(x['hh_size'], x['hh_frac']) for x in self.data] \
            if (st is None or et is None) \
            else [(x['hh_size'], x['hh_frac']) for x in self.data.where('(t>=st) & (t<et)')]
        for cur_case in size_hh_frac:
            # store ages for hh of size (1..6, 7+)
            hh_fracs[min(cur_case[0], 7)].append(cur_case[1])
        avg_hh_fracs = defaultdict(list)
        for cur_hh_size in range(1, 8):
            avg_hh_fracs[cur_hh_size].append(
                NaN if cur_hh_size not in hh_fracs
                else np.mean(hh_fracs[cur_hh_size]))
        return avg_hh_fracs

    def get_hh_fracs_by_age(self, min_age, max_age, t_step):
        """
        Not currently used...

        :param min_age:
        :param max_age:
        :param t_step:
        :return:
        """
        hh_fracs_by_age = defaultdict(list)
        for x in self.data:
            if min_age <= x['age'] < max_age:
                hh_fracs_by_age[x['t'] // t_step].append(x['hh_frac'])
        return sorted([(k, np.mean(v), np.std(v)) for k, v in list(hh_fracs_by_age.items())])

    #########################################################################
    ## clustering of susceptibility in households

    def get_hh_at_risk_clustering(self, st=None, et=None):
        """
        Get the mean number of susceptibles in a household at the point of disease introduction
        by household size over time.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :return: A list of tuples, containing mean/SD age of first infection for that time period.
        """

        # set up object for storing by hh_size
        hh_at_risk_by_size = [[] for _ in range(8)]

        # process introductory cases occurring in this time period
        for x in self.data.where('(t>=st) & (t<et)'):
            if x['hh_at_risk'] > 0:
                hh_at_risk_by_size[min(x['hh_size'], 7)].append(int(x['hh_at_risk']))

        for i, x in enumerate(hh_at_risk_by_size):
            if len(x) == 0:
                hh_at_risk_by_size[i].append(0)

        return [np.mean(x) for x in hh_at_risk_by_size if len(x) > 0]

    def get_hh_at_risk_prop(self, st=None, et=None):
        """
        Get the proportion of households containing more than one susceptible at the time of
        disease introduction by household size over time.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :return: A numpy array of susceptible proportions by household size
        """

        # set up object for storing by hh_size
        hh_at_risk_total = np.zeros(8)
        hh_at_risk_cluster = np.zeros(8)  # more than one susceptible in hh

        # process introductory cases occurring in this time period
        for x in self.data.where('(t>=cur_st) & (t<cur_et)'):
            if x['hh_at_risk'] > 0:
                hh_at_risk_total[min(x['hh_size'], 7)] += 1
                if x['hh_at_risk'] > 1:
                    hh_at_risk_cluster[min(x['hh_size'], 7)] += 1

        return hh_at_risk_cluster / hh_at_risk_total

    #########################################################################
    ## age distributions

    def get_case_age_dist(self, st=None, et=None, state=None, age_bins=None, first=False):
        """
        Get a distribution of the ages (in years) at which cases occurred.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :param state: The disease state to get distribution over
        :type state: string
        :param age_bins: The age bins to aggregate over (in years).
        :type age_bins: list
        :param first: If `True`, restrict to first cases.
        :type first: bool
        :return: A list containing the distribution of case counts over age classes
        """

        if not state:
            state = self.entry_state

        if not age_bins:
            age_bins = list(range(101))

        if st is None or et is None:
            ages = [x['age'] for x in self.data if (x['state'] == state) and (x['first'] if first else True)]
        else:
            ages = [x['age'] for x in self.data.where('(t>=st) & (t<et)') if (x['state'] == state) and
                                                                             (x['first'] if first else True)]
        age_dist = np.histogram(ages, bins=age_bins)[0]
        return age_dist

    def get_case_age_dist_days(self, st=None, et=None, state=None, age_days=None, max_age=2, first=False):
        """
        Get a distribution of the ages (in days) at which cases occurred.

        Default bins are the first 12 months.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :param state: The disease state to get distribution over
        :type state: string
        :param age_days: The age bins to aggregate over (in days).
        :type age_days: list
        :param max_age: The maximum age to aggregate up to
        :type max_age: int
        :param first: If `True`, restrict to first cases.
        :type first: bool
        :return: A list containing the distribution of case counts over age classes
        """

        if not state:
            state = self.entry_state

#        if not age_days:
#            age_days = month_days
#        age_bins = []
#        for cur_age in range(max_age):
#            age_bins.extend([(cur_age) * 364 + x for x in age_days])
#        age_bins.extend([364*max_age])
#
#        max_age = (age_bins[-1] / 364) + 1

        age_bins = age_days
        print(age_bins)

        if st is None or et is None:
            ages = [x['age'] * 364 + x['age_days']
                    for x in self.data
                    if (x['age'] < max_age) and (x['state'] == state) and (x['first'] if first else True)]
        else:
            #ages = [x['age_days'] for x in self.data.where('(t>=st) & (t<et)') if x['age'] < max_age]
            ages = [(x['age'] * 364 + x['age_days'])
                    for x in self.data.where('(t>=st) & (t<et)')
                    if (x['age'] < max_age) and (x['state'] == state) and (x['first'] if first else True)]

        age_dist = np.histogram(ages, bins=age_bins)[0]
        print(age_dist)
        return age_dist

    #########################################################################
    ## hh source proportion by age over time

    def get_hh_infections_days(self, st=None, et=None, state=None, age_days=None, max_age=None):
        """
        Get proportion of infant infections due to household sources, by age.
        """

        if not state:
            state = self.entry_state

#        if not age_days:
#            age_days = month_days
#        age_bins = []
#        for cur_age in range(max_age):
#            age_bins.extend([(cur_age) * 364 + x for x in age_days])
#        age_bins.extend([364*max_age])
#
#        max_age = (age_bins[-1] / 364) + 1

        age_bins = age_days

        if st is None or et is None:
            ages = [x['age'] * 364 + x['age_days']
                    for x in self.data
                    if (x['age'] < max_age) and (x['state'] == state) and (x['hh_source'])]
        else:
            ages = [x['age'] * 364 + x['age_days']
                    for x in self.data.where('(t>st) & (t<et)')
                    if (x['age'] < max_age) and (x['state'] == state) and (x['hh_source'])]

        age_dist = np.histogram(ages, bins=age_bins)[0]
        return age_dist

    def get_hh_source_age_days(self, st, et, state, age_days, max_age=None):

#        age_bins = []
#        for cur_age in range(max_age):
#            age_bins.extend([(cur_age) * 364 + x for x in age_days])
#        age_bins.extend([364*max_age])

        age_bins = age_days
        print(age_bins)
        #max_age = (age_bins[-1] / 364) + 1

        avg_ages = []
        for cur_min_age, cur_max_age in zip(age_bins[:-1], age_bins[1:]):
            ages = [x['source_age']
                    for x in self.data.where('(t>st) & (t<et)')
                    if (cur_min_age <= (x['age'] * 364 + x['age_days']) < cur_max_age) and (x['state'] == state) and (x['hh_source'])]
            avg_ages.append(np.mean(ages))
        return avg_ages

    def bin_age_dist(self, age_dist, age_bins):
        """
        Convert an distribution of age counts into a distribution of binned counts.
        Can be used for both years and months/days

        :param age_dist: The age distribution to bin
        :type age_dist: list
        :param age_bins: The bin edges
        :type age_bins: list
        :return: A binned age distribution (list)
        """

        binned_dist = []
        for cur_bin_l, cur_bin_r in zip(age_bins, age_bins[1:] + [101]):
            binned_dist.append(sum(age_dist[cur_bin_l:cur_bin_r]))
        return binned_dist

    def get_mean_age_dist_bin(self, age_bins, st, et):
        """
        Get mean binned age distribution over specified period.

        :param age_bins: The bin edges
        :type age_bins: list
        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        """

        binned_age_dists = []
        for x in self.age_dists.where('(t>=st) & (t<et)'):
            binned_age_dists.append(self.bin_age_dist(x['dist'], age_bins))
        return np.mean(binned_age_dists, 0)

    def get_mean_age_dist(self, st, et):
        """
        Get the mean age distribution over the specified period.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :return: The binned age distribution.
        """

        return np.mean([x['dist'] for x in self.age_dists.where('(t>=st) & (t<et)')], 0)

    def get_mean_age_dist_month(self, st, et):
        """
        Get the mean (month) household size distribution over the specified period.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :return: The binned age distribution.
        """

        return np.mean([x['dist'] for x in self.age_dists_month.where('(t>=st) & (t<et)')], 0)

    def get_mean_age_dist_month_bin(self, age_bins_month, st, et):
        """
        Get mean binned monthly age distribution over specified period.

        :param age_bins: The bin edges (NB: in months, not days)
        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :return: The binned age distribution.
        """

        binned_age_dists = []
        for x in self.age_dists_month.where('(t>=st) & (t<et)'):
            binned_age_dists.append(self.bin_age_dist(x['dist'], age_bins_month))
        return np.mean(binned_age_dists, 0)

    #########################################################################
    ## household size distributions

    def get_mean_hh_size_dist(self, st, et):
        """

        Get the mean age distribution over the specified period.

        :param st: The start time step
        :type st: int
        :param et: The end time step
        :type et: int
        :return: The binned age distribution.
        """
        return np.mean([x['dist'] for x in self.hh_size_dists.where('(t>=st) & (t<et)')], 0)

    #########################################################################
    ## overall incidence

    def get_annual_incidence(self, times, dt, loc=None):
        """
        Get incidence binned over requested time periods, potentially in particular location type.

        :param times:
        :param dt:
        :param loc:
        :return:
        """
        sizes = self.get_mean_pop_size_time(times['st'], times['et'], dt)
        cases = self.get_new_case_counts_by_time(times['st'], times['et'], dt, loc)
        return np.array(cases) / np.array(sizes) #/ dt * times['t_per_year']  # last bit makes it annual (?)

    def get_first_infection_incidence(self, times, dt):
        """
        Get incidence of first infections binned over requested time periods.

        :param times:
        :param dt:
        :return:
        """
        sizes = self.get_mean_pop_size_time(times['st'], times['et'], dt)
        cases = self.get_new_case_counts_by_time_first(times['st'], times['et'], dt)
        return np.array(cases) / np.array(sizes) #/ dt * times['t_per_year']  # last bit makes it annual (?)


    #########################################################################
    ## incidence by age

    def get_incidence_by_age_infant(self, age_dist, st, et, max_age):
        """
        Infant version of below; note curently doesn't handle multiple infectious states.
        Could be extended to do so by adding further parameter to get_case_age_dist...
        """

        filtered = self.filter_data_by_state()
        incs = {}
        for k, v in list(filtered.items()):
            ages = self.get_case_age_dist_days(st, et, k, max_age=max_age)
            duration = max((et - st) / float(self.t_per_year), 1.0)
            incs[k] = [(float(x) / (y * duration) if y > 0 else 0.0) for x, y in zip(ages, age_dist)]
        return incs


    def get_hh_infections_infant(self, st, et, max_age):

        filtered = self.filter_data_by_state()
        props = {}
        for k, v in list(filtered.items()):
            ages = self.get_case_age_dist_days(st, et, k, max_age=max_age)
            hh_ages = self.get_hh_infections_days(st, et, k, max_age=max_age)
            props[k] = [(float(x) / y if y > 0 else 0.0) for x, y in zip(hh_ages, ages)]
        return props



#    def get_incidence_by_age_month(self, age_dist, st, et, max_age, first=False):
#        """
#        Get incidence by age, broken down by months up to max_age.
#        """
#
#        filtered = self.filter_data_by_state()
#        incs = {}
#        for k, v in filtered.items():
#            ages = self.get_case_age_dist_days(st, et, k, max_age=max_age, first=first)
#            # NB: following weirdness is to deal with fact that age_dists were stored at
#            # coarser than desired resolution -- should be able to remove.
#            print "len ages:", len(ages)
#            print "len age_dist:", len(age_dist)
#            factor = len(ages) / len(age_dist)
#            ext_age_dist = [x for x in chain(*[[x/float(factor)]*factor for x in age_dist])]
#            print "len ext_age_dist:", len(ext_age_dist)
#            print ext_age_dist[:10]
#            duration = max((et - st) / float(self.t_per_year), 1.0)
#            incs[k] = [(float(x) / (y * duration) if y > 0 else 0.0) for x, y in zip(ages, ext_age_dist)]
#        return incs

    def get_incidence_by_age(self, age_dist, st, et, first=False):
        """
        Get annual incidence broken down by age, calculated over the 
        specified interval.  Age distribution should be specified in 
        non-normalised form (i.e., number of individuals in each class,
        rather than proportion).

        Returns incidences in a dictionary keyed by actual disease state, to
        allow for handling of multiple infectious states (e.g., primary and
        secondary).

        TODO: 
        * probably not flexible enough to handle multi-year bins atm (?)
        * something else...
        """

        filtered = self.filter_data_by_state()
        incs = {}
        for k, v in list(filtered.items()):
            ages = self.get_case_age_dist(st, et, k, first=first)
            duration = max((et - st) / float(self.t_per_year), 1.0)
            incs[k] = [(float(x) / (y * duration) if y > 0 else 0.0) for x, y in zip(ages, age_dist[:-1])]
        return incs

    def get_incidence_by_age_snapshots(self, st, et, num_snapshots, first=False):
        """
        Get series of incidence by age distributions at evenly spaced time intervals
        between st and et.
        """
        dt = (et - st) // num_snapshots
        inc_age_dists = defaultdict(list)
#        print num_snapshots, dt
        for i in range(num_snapshots):
            cur_st = st + (i * dt)
            cur_et = cur_st + dt
            age_dist = self.get_mean_age_dist(cur_st, cur_et)
            cur_inc = self.get_incidence_by_age(age_dist, cur_st, cur_et, first)
            for k, v in list(cur_inc.items()):
                inc_age_dists[k].append(v)
        return inc_age_dists

#    def get_incidence_by_age_snapshots_month(self, st, et, num_snapshots, first=False):
#        dt = (et - st) / num_snapshots
#        inc_age_dists = []
##        print "num snapshots:", num_snapshots, "; dt:", dt
#        for i in range(num_snapshots):
#            cur_st = st + (i * dt)
#            cur_et = cur_st + dt
#            age_dist = self.get_mean_age_dist_month(cur_st, cur_et)
#            cur_inc = self.get_incidence_by_age_infant_month(age_dist, cur_st, cur_et, 5, first)
#            inc_age_dists.append(cur_inc[self.entry_state])
#        return inc_age_dists

    def get_incidence_by_age_infant_snapshots(self, st, et, num_snapshots):
        """
        NB: this seems a bit broken!!! (fixed now)
        """
        dt = (et - st) // num_snapshots
        inc_age_infant_dists = []
        for i in range(num_snapshots):
            cur_st = st + (i * dt)
            cur_et = cur_st + dt
            age_dist_quarters = self.get_mean_age_dist_month(st, et)
            print(("A", age_dist_quarters))
            cur_inc = self.get_incidence_by_age_infant(age_dist_quarters, cur_st, cur_et, 2)
            inc_age_infant_dists.append(cur_inc[self.entry_state])
        return inc_age_infant_dists

    def get_incidence_by_age_bin(self, age_bins, st, et, dt):
        """
        Get incidence by age curves, for specified age bins, over specified intervals
        (st, et = start time, end time; dt = time step size).
        
        Uses closest available stored age distribution.
        """
        incs = defaultdict(list)
        print((st, et, dt))
        duration = float(dt) / self.t_per_year
        for cur_st in range(st, et, dt):
            cur_et = cur_st + dt
            age_dist = self.get_mean_age_dist_bin(age_bins, cur_st, cur_et)
            filtered = self.filter_data_by_state(cur_st, cur_et)
            for k, v in list(filtered.items()):
                ages = []
                for r in v:
                    ages.append(r['age'])
                d, b = np.histogram(ages, bins=age_bins + [101], normed=False)
                incs[k].append([(float(x) / (y * duration)) if y > 0 else 0 for x, y in zip(d, age_dist)])
        return incs

    def get_incidence_by_age_bin_months(self, age_bins_month, age_days, st, et, dt, max_age=5):
        """
        Get incidence by age curves, for specified age bins, over specified intervals
        (st, et = start time, end time; dt = time step size).

        Uses closest available stored age distribution.
        """
        incs = defaultdict(list)
        duration = float(dt) / self.t_per_year
        for cur_st in range(st, et, dt):
            cur_et = cur_st + dt
            age_dist = self.get_mean_age_dist_month_bin(age_bins_month, cur_st, cur_et)
            filtered = self.filter_data_by_state(cur_st, cur_et)
            for k, v in list(filtered.items()):
                ages = self.get_case_age_dist_days(cur_st, cur_et, k, age_days=age_days, max_age=max_age)
                #ages = []
                #for r in v:
                #    ages.append(r['age'])
                #d, b = np.histogram(ages, bins=age_bins, normed=False)
                incs[k].append([(float(x) / (y * duration)) if y > 0 else 0 for x, y in zip(ages, age_dist)])
        return incs

    def get_hh_infections_by_age_bin_months(self, age_bins, st, et, dt, max_age=2):

        props = defaultdict(list)
        for cur_st in range(st, et, dt):
            cur_et = cur_st + dt
            filtered = self.filter_data_by_state(cur_st, cur_et)
            for k, v in list(filtered.items()):
                ages = self.get_case_age_dist_days(cur_st, cur_et, k, age_days=age_bins, max_age=max_age)
                hh_ages = self.get_hh_infections_days(cur_st, cur_et, k, age_days=age_bins, max_age=max_age)
                props[k].append([(float(x) / y if y > 0 else 0.0) for x, y in zip(hh_ages, ages)])
        return props

    def get_hh_source_age_by_age_bin_months(self, age_bins, st, et, dt, max_age=2):

        avg_ages = defaultdict(list)
        for cur_st in range(st, et, dt):
            cur_et = cur_st + dt
            filtered = self.filter_data_by_state(cur_st, cur_et)
            for k, v in list(filtered.items()):
#                ages = self.get_case_age_dist_days(cur_st, cur_et, k, age_days=age_bins, max_age=max_age)
                hh_source_ages = self.get_hh_source_age_days(cur_st, cur_et, k, age_days=age_bins, max_age=max_age)
                avg_ages[k].append(hh_source_ages)
        return avg_ages

    #########################################################################
    ## incidence by household size (and type)

    def get_case_hh_size_distribution(self, max_size, st=None, et=None, normed=False):
        """
        Returns the distribution of household sizes for all observed cases.

        Used in FF100 sims?
        """
        if st is not None:
            hh_sizes = [x['hh_size'] for x in self.data.where('(t>=st) & (t<et)')]
        else:
            hh_sizes = [x['hh_size'] for x in self.data]
        # bins: specify the right edge of the final bin; starts at 1 (no hh of size 0)
        d, b = np.histogram(hh_sizes, bins=np.arange(0, max_size + 1) + 1, normed=normed)
        return d

    def get_case_hh_size_dist_hh(self, max_size, num_cases=None, st=None, et=None, normed=False):
        """
        Returns the distribution of household sizes for households with at least one observed case.
        """
        hhs_by_size = [set() for _ in range(max_size+1)]
        if st is not None:
            for cur_case in self.data[:num_cases].where('(t>=st) & (t<et)'):
                hhs_by_size[min(cur_case['hh_size'], max_size)].add(cur_case['hh_id'])
        else:
            for cur_case in self.data[:num_cases]:
                hhs_by_size[min(cur_case['hh_size'], max_size)].add(cur_case['hh_id'])

        counts = np.array([len(hhs_by_size[x]) for x in range(1, max_size+1)], dtype=np.float)
        if normed:
            counts /= np.sum(counts)
        return counts

    def get_incidence_by_hh_size(self, ind_hh_dist, st, et, max_size):
        """
        Get annual incidence broken down by household size.  ind_hh_dist is 
        distribution of individuals by household size (again, non-normalised).
        
        max_size: aggregate all households of this size or above
        """
        filtered = self.filter_data_by_state()
        incs = {}
        for k, v in list(filtered.items()):
            hh_sizes = [x['hh_size'] for x in v if st < x['t'] < et]
            d, b = np.histogram(hh_sizes, bins=len(ind_hh_dist),
                                range=(0, len(ind_hh_dist)), normed=False)
            binned_d = list(d[0:max_size]) + [sum(d[max_size:])]
            binned_dist = list(ind_hh_dist[:(max_size - 1)]) + [sum(ind_hh_dist[max_size - 1:])]
            duration = max((et - st) / float(self.t_per_year), 1.0)
            #incs[k] = [(float(x) / ((y[0] + 1) * y[1] * duration) if y[1] > 0 else 0.0) for x, y in
            #           zip(binned_d, enumerate(binned_dist))]
            incs[k] = [(float(x) / (y[1] * duration) if y[1] > 0 else 0.0) for x, y in
                       zip(binned_d, enumerate(binned_dist))]
        return incs

    def get_incidence_by_hh_size_snapshots(self, st, et, num_snapshots, max_hh_size=10):
        """
        NB: this seems a bit broken!!! (fixed now)
        """
        dt = (et - st) / num_snapshots
        inc_hh_dists = []
        for i in range(num_snapshots):
            cur_st = st + (i * dt)
            cur_et = cur_st + dt
            hh_size_dist = self.get_mean_hh_size_dist(cur_st, cur_et)
            cur_inc = self.get_incidence_by_hh_size(hh_size_dist, cur_st, cur_et, max_hh_size)
            inc_hh_dists.append(cur_inc[self.entry_state])
        return inc_hh_dists

    def get_incidence_by_hh_size_time(self, st, et, dt):
        incs = []
        for cur_st in range(st, et, dt):
            cur_et = cur_st + dt
            hh_size_dist = self.get_mean_hh_size_dist(cur_st, cur_et)
            cur_incs = self.get_incidence_by_hh_size(hh_size_dist, cur_st, cur_et, 10)
            incs.append(cur_incs[self.entry_state])
        return incs

    def get_incidence_by_hh_type(self, ind_hh_type_dist, st, et):
        """
        Get annual incidence broken down by household type.  ind_hh_type_dist 
        is distribution of individuals by household type 
        (again, non-normalised).
        """

        filtered = self.filter_data_by_state()
        incs = {}
        for k, v in list(filtered.items()):
            hh_types = ['kids' if x['hh_type'].endswith('kids') else 'nokids'
                        for x in v if st < x['t'] < et and 18 < x['age']]
            counts = (hh_types.count('kids'), hh_types.count('nokids'))
            duration = max((et - st) / float(self.t_per_year), 1.0)
            incs[k] = [float(x) / (y * duration) for x, y in \
                       zip(counts, ind_hh_type_dist)]
        return incs

    ###############################################
    # secondary household attach proportion

    def get_shap(self, st=None, et=None):
        """
        NB: SHAP calculation should only be done for index cases in households (otherwise will be underestimating...)
        """
        shaps_true = []
        shaps_obs = []
        shaps_true_first = []
        shaps_obs_first = []
        data = self.data if st is None else self.data.where('(t>=st) & (t<et)')
        for x in data:
            #print x, x['hh_infections'], x['hh_at_risk'] - 1
            if x['hh_size'] < 2: continue
            if not x['hh_index_case']: continue
            shaps_true.append(float(x['hh_infections']) / (x['hh_at_risk'] - 1) if x['hh_at_risk'] > 1 else 0)
            shaps_obs.append(float(x['hh_infections']) / (x['hh_size'] - 1) if x['hh_size'] > 1 else 0)
            if not x['hh_index_case_first']: continue
            shaps_true_first.append(float(x['hh_infections']) / (x['hh_at_risk'] - 1) if x['hh_at_risk'] > 1 else 0)
            shaps_obs_first.append(float(x['hh_infections']) / (x['hh_size'] - 1) if x['hh_size'] > 1 else 0)

        return np.mean(shaps_true), np.mean(shaps_obs), np.mean(shaps_true_first), np.mean(shaps_obs_first)

    ################################################

    def estimate_R0(self):
        case_t = self.get_new_data_by_time()

        cases = list(case_t.items())
        cases.sort()
        xvals, yvals = list(zip(*cases))

        xvals = np.array(xvals) - 36400
        yvals = np.array(yvals)

        func = lambda x, b, c: np.exp(b * x) + c
        popt, pcov = curve_fit(func, xvals, yvals)

        for i in range(2, len(xvals) + 1):
            popt, pcov = curve_fit(func, xvals[:i], yvals[:i])
            slope = popt[0]
            #print i, ':', slope

    #########################################################################
    ## demographic stuff (that probably ought be in its own observer)

    def get_pop_size_time(self):
        time_values = [x['t'] for x in self.age_dists]
        pop_sizes = [sum(x['dist']) for x in self.age_dists]
        return time_values, pop_sizes

    def get_mean_pop_size_time(self, st, et, dt):
        sizes = []
        for cur_st in range(st, et, dt):
            cur_et = cur_st + max(dt, self.t_per_year / 4)
            pop_sizes = [sum(x['dist']) for x in self.age_dists.where('(t >= cur_st) & (t < cur_et)')]
            sizes.append(np.mean(pop_sizes))
        return sizes

    def get_age_dists_time(self, t_snapshots=None):
        age_dists = []
        if not t_snapshots:
            t_snapshots = self.age_dists.col('t')

        for cur_t in t_snapshots:
            raw_dist = [x['dist'] for x in self.age_dists.where('t==cur_t')][0]
            total = sum(raw_dist)
            prop_dist = [x / float(total) for x in raw_dist]
            age_dists.append(prop_dist)
        return age_dists

    def get_hh_size_dists_time(self, t_snapshots=None, max_hh_size=10):
        hh_size_dists = []
        if not t_snapshots:
            t_snapshots = self.age_dists.col('t')

        for cur_t in t_snapshots:
            raw_dist = [x['dist'] for x in self.hh_size_dists.where('t==cur_t')][0]
            total = sum(raw_dist)
            prop_dist = [x / float(total) for x in raw_dist[:(max_hh_size - 1)]] + [
                sum(raw_dist[(max_hh_size - 1):]) / float(total)]
            hh_size_dists.append(prop_dist)
        return hh_size_dists

    #########################################################################
    ## utility functions

    def get_snapshot_times(self, num_snapshots):
        time_values = [x['t'] for x in self.age_dists]
        length = float(len(time_values))
        t_snapshots = [time_values[int(ceil(x * length / num_snapshots))] for x in range(num_snapshots)]
        return t_snapshots

# mydata = np.array([(106869, 7, 4, 5, np.array([201032, 201813, 169339, 204419, 175358]))],
#               dtype=[('hh.id', 'i'), ('N', 'i'), ('det.num', 'i'), ('fs', 'i'), ('id', 'object')])

    def build_hh_structure(self):
        hh_count = len(set(self.data.col('hh_id')))
        print(hh_count)
        hh_data = {}  # maps hh_id to hh_struct dict

        second_intro = {}

        for x in self.data:
            cur_hh_id = x['hh_id']

            if cur_hh_id in second_intro:
                continue

            if x['hh_index_case'] and not x['hh_index_case_first']:
                print((x['hh_id']))
                second_intro[cur_hh_id] = True
                continue

            if cur_hh_id not in hh_data:
                hh_data[cur_hh_id] = {
                    'hh_id': x['hh_id'],
                    'N': x['hh_size'],
                    'det_num': 0,
                    'fs': [],
                    'id': [],
                    'inf_id': [],
                    'all_times': [],
                    'det_times': [],
                    'ext_flag': [],
                    'ser_fs': 0,
                    'W': []
                }
            hh_data[cur_hh_id]['fs'] = 0
            hh_data[cur_hh_id]['inf_id'].append(x['source_id'])
            hh_data[cur_hh_id]['all_times'].append(x['t'])
            hh_data[cur_hh_id]['ext_flag'].append(x['hh_source'])
        return hh_data




    def export_as_csv(self, filename):
        subset = [self.data.col('t'), self.data.col('t_E'), self.data.col('ID'), self.data.col('age'),
                  self.data.col('hh_id'), self.data.col('hh_size'), self.data.col('hh_source'),
                  self.data.col('source_id')] #, self.data.col('ext_contacts'), self.data.col('ext_infections')]
        s_array = np.array(subset)
        s_array = np.transpose(s_array)
        np.savetxt(filename, s_array, delimiter=',', fmt='%d')

    def export_episodes_as_csv(self, filename):
        subset = [self.data.col('t') / self.t_per_year, (self.data.col('t_R') - self.data.col('t')) / self.t_per_year,
                  self.data.col('ID'), self.data.col('hh_id'), self.data.col('age'), self.data.col('hh_size')]
        s_array = np.array(subset)
        s_array = np.transpose(s_array)
        np.savetxt(filename, s_array, delimiter=',')

    def write_age_incidence(self, filename, times, first):

        st = times['et'] - (10*self.t_per_year)
        age_inc = self.get_incidence_by_age_snapshots(
            st, times['et'], (times['et']-st) / self.t_per_year, first=first)
        np.savetxt(filename, np.array(age_inc[self.entry_state]), delimiter=',')

#    def write_age_incidence_infant_month(self, filename, times, first):
#
#        st = times['et'] - (10*self.t_per_year)
#        age_inc = self.get_incidence_by_age_snapshots_month(
#            st, times['et'], (times['et']-st) / self.t_per_year, first=first)
#        np.savetxt(filename, np.array(age_inc), delimiter=',')

    #########################################################################
    ### OUTPUT functions ###
    def output_all(self, p, times):

        return

        # - #  -- demography outputs -- # - # - # - # - # - # - # - # - # - # - #

        # A hack to dump the population size (obviously only of interest in a growing population!)
        # eventually this should be handled by a demography observer that would also store household
        # and age distributions (and could potentially be used by case observer for denominator values)

        x_offset = 1910-(times['st'] / self.t_per_year)

        #if p['dyn_rates']:
        print("Generating demographic plots")
        output_timeseries(os.path.join(p['prefix'], 'demog_pop_size.png'), times,
                          [self.get_pop_size_time()[1]], 'Population Size',
                          x_offset=x_offset)

        t_snapshots = self.get_snapshot_times(2)
        age_dists = self.get_age_dists_time(t_snapshots)

        num = 0.0
        denom = 0.0
        for cur_age, cur_count in enumerate(age_dists[0]):
            num += cur_age * cur_count
            denom += cur_count
        print((num / denom))

#        plot_statistics(os.path.join(p['prefix'], 'demog_age_dists.png'), range(len(age_dists[0])),
#                        age_dists, 'Age', 'Proportion',
#                        series_labels=['t=%d' % (x / p['t_per_year']) for x in t_snapshots])
#
#        hh_size_dists = self.get_hh_size_dists_time(t_snapshots, max_hh_size=10)
#        plot_statistics(os.path.join(p['prefix'], 'demog_hh_size_dists.png'), range(1, len(hh_size_dists[0]) + 1),
#                        hh_size_dists, 'Household size', 'Proportion',
#                        series_labels=['t=%d' % (x / p['t_per_year']) for x in t_snapshots])
#
#        age_bins = [0, 1, 5, 10, 20, 40, 60]
#        age_dists_all = zip(*[self.bin_age_dist(x, age_bins) for x in self.get_age_dists_time()])
#
#        output_stacked_bars_timeseries(os.path.join(p['prefix'], 'demog_age_props.png'),
#                                       times, age_dists_all, 'Proportion of population',
#                                       labels=['<%d' % age_bins[0]] + ['%d-%d' % (x, y) for x, y in
#                                                                       zip(age_bins[1:], age_bins[2:] + [100])])
#
#        output_stacked_bars_timeseries(os.path.join(p['prefix'], 'demog_hh_size_props.png'),
#                                       times, zip(*self.get_hh_size_dists_time(max_hh_size=7)), 'Proportion',
#                                       labels=['%d' % x for x in range(1, 7)] + ['7+'])

        # - #  -- timeseries outputs -- # - # - # - # - # - # - # - # - # - # - #


#        # output case count per time period
#        # TODO: make options for daily, weekly, etc.
#        print "Generating case timeseries plot"
#        output_timeseries(os.path.join(p['prefix'], 'cases_new.png'), times,
#                          self.get_new_case_counts_by_time(
#                              times['st'], times['et'], 1),
#                          'New cases', x_offset=x_offset)
#
#        output_timeseries(os.path.join(p['prefix'], 'cases_cum.png'), times,
#                          self.get_cum_case_count(times['st'], times['et'], self.t_per_year),
#                          'Cumulative cases', x_offset=x_offset)
#
#
#        incs_total = self.get_annual_incidence(times, 1)
#        incs_first = self.get_first_infection_incidence(times, 1)
#
#        output_timeseries(os.path.join(p['prefix'], 'cases_weekly_incidence.png'), times,
#                          [incs_total, incs_first], 'Weekly Incidence', labels=['All infections', 'First infections'],
#                          x_offset=x_offset)
#
##        output_timeseries(os.path.join(p['prefix'], 'cases_weekly_incidence_ratio.png'), times,
##                          np.array(incs_first) / np.array(incs_total), 'Ratio of first to total infections',
##                          x_offset=x_offset)
#
#        if (times['et'] - times['st']) / self.t_per_year > 10:
#            print "Generating annual incidence plot"
#            incs_total = self.get_annual_incidence(times, self.t_per_year)
#            incs_first = self.get_first_infection_incidence(times, self.t_per_year)
#            output_timeseries(os.path.join(p['prefix'], 'cases_total_incidence.png'), times,
#                              [incs_total, incs_first], 'Annual Incidence',
#                              labels=['All infections', 'First infections'], x_offset=x_offset)
##
##        #self.write_age_incidence(os.path.join(p['prefix'], 'age_incidence_all.csv'), times, first=False)
##        #self.write_age_incidence(os.path.join(p['prefix'], 'age_incidence_first.csv'), times, first=True)
##
##
##
##        print "Writing age incidence infant CSV output"
##        self.write_age_incidence_infant_month(os.path.join(p['prefix'], 'age_incidence_first_infant.csv'),
##                                              times, first=True)
#
#        # - #  -- incidence outputs -- # - # - # - # - # - # - # - # - # - # - #
#
#        # age incidence
#        print "Generating incidence by age plots"
#        #age_dist = self.get_mean_age_dist(times['st'], times['et'])
#        #age_inc_dat = self.get_incidence_by_age(age_dist, times['st'], times['et'])
#
#        #incs_snapshots = [age_inc_dat[self.entry_state]]
#        #incs_snapshots = []
#
#        #if p['dyn_rates']:
#        incs_snapshots = self.get_incidence_by_age_snapshots(times['st'], times['et'], 10)
#
#        print "number of snapshots: ", len(incs_snapshots['In'])
#        for k, v in incs_snapshots.items():
#            plot_statistics(os.path.join(p['prefix'], 'cases_incidence_by_age_snapshots_%s.png' % k),
#                            range(1, len(v[0])+1), v, 'Age', 'Incidence', xlogscale=True,
#                            series_labels=['t=%d' % (x_offset + (x / p['t_per_year'])) for x in t_snapshots],
#                            legend_loc='upper right', legend_title='Time period beginning:')
#
#        # age incidence over time (years)
#
#        age_bins = [0, 1, 5, 10, 20, 40, 60]
#        incs = self.get_incidence_by_age_bin(age_bins, times['st'], times['et'], times['t_per_year'])
#
#        for k, v in incs.items():
#            output_timeseries(os.path.join(p['prefix'], 'cases_incidence_by_age_time_%s.png' % k),
#                              times, np.transpose(np.array(v)), 'Incidence',
#                              labels=['%d-%d' % (x, y) for x, y in zip(age_bins, age_bins[1:] + [100])],
#                              x_offset=x_offset, legend_title='Age group (years)')
#
#        output_timeseries(os.path.join(p['prefix'], 'cases_incidence_by_age_time_total.png'),
#                          times, np.transpose(np.sum(np.array([v for v in incs.values()]), axis=0)), 'Incidence',
#                          labels=['%d-%d' % (x, y) for x, y in zip(age_bins, age_bins[1:] + [100])],
#                          x_offset=x_offset, legend_title='Age group (years)')
#
        # age incidence over time (months)

#        age_bins_months = [0, 2, 4, 6, 12, 18, 24, 60]
#        age_bins_days = [(x/12) * 364 + month_days[x%12] for x in age_bins_months]
#        incs = self.get_incidence_by_age_bin_months(age_bins_months, age_bins_days, times['st'], times['et'], times['t_per_year'])
#
#        for k, v in incs.items():
#            output_timeseries(os.path.join(p['prefix'], 'cases_incidence_by_age_time_months_%s.png' % k),
#                              times, np.transpose(np.array(v)), 'Incidence',
#                              labels=['%d-%d' % (x, y) for x, y in zip(age_bins_months[:-1], age_bins_months[1:])],
#                              x_offset=x_offset, legend_title='Age group (months)')
#
#        output_timeseries(os.path.join(p['prefix'], 'cases_incidence_by_age_time_months_total.png'),
#                          times, np.transpose(np.sum(np.array([v for v in incs.values()]), axis=0)), 'Incidence',
#                          labels=['%d-%d' % (x, y) for x, y in zip(age_bins_months, age_bins_months[1:] + [60])],
#                          x_offset=x_offset, legend_title='Age group (months)')
#
#        age_bins_months = [0, 364, 728, 364*5]
#        props = self.get_hh_infections_by_age_bin_months(age_bins_days, times['st'], times['et'],
#                                                         times['t_per_year']*10, max_age=5)
#        for k, v in props.items():
#            output_timeseries(os.path.join(p['prefix'], 'hh_source_props_by_age_time_months_%s.png' % k),
#                              times, np.transpose(np.array(v)), 'Proportion',
#                              labels=['%d-%d' % (x, y) for x, y in zip(age_bins_months[:-1], age_bins_months[1:])],
#                              x_offset=x_offset, legend_title='Age group (months)')
#
        age_bins_months = [0, 364, 728, 364*5]
        avg_source_ages = self.get_hh_source_age_by_age_bin_months(age_bins_months, times['st'], times['et'],
                                                                   times['t_per_year']*10, max_age=5)

        for k, v in list(avg_source_ages.items()):
            print(v)
            print((len(v)))
            output_timeseries(os.path.join(p['prefix'], 'hh_source_age_by_age_time_months_%s.png' % k),
                              times, np.transpose(np.array(v)), 'Average age of source (years)',
                              labels=['%d-%d' % (x, y) for x, y in zip(age_bins_months, age_bins_months[1:] + [364*5])],
                              x_offset=x_offset, legend_title='Age group (months)')




#        return

        # household size incidence
        print("Generating incidence by hh size plots")
        max_hh_size = 7
        #hh_size_dist = self.get_mean_hh_size_dist(times['st'], times['et'])
        #hh_inc_dat = self.get_incidence_by_hh_size(hh_size_dist, times['st'], times['et'], max_hh_size + 1)

        # start with overall incidence
        #incs_snapshots = [hh_inc_dat[self.entry_state]]

        #if p['dyn_rates']:
        # add on incidence across evenly spaced intervals
        incs_snapshots = self.get_incidence_by_hh_size_snapshots(times['st'], times['et'], 10,
                                                                 max_hh_size=max_hh_size + 1)

        plot_statistics(os.path.join(p['prefix'], 'cases_incidence_by_hh_size_snapshots.png'),
                        list(range(len(incs_snapshots[0]))), incs_snapshots, 'Household size', 'Incidence',
                        series_labels=['t=%d' % (x_offset + (x / p['t_per_year'])) for x in t_snapshots],
                        legend_loc='upper left', legend_title='Time period beginning:')

        incs = self.get_incidence_by_hh_size_time(times['st'], times['et'], times['t_per_year'])

        output_timeseries(os.path.join(p['prefix'], 'cases_incidence_by_hh_size_time.png'),
                          times, np.transpose(np.array(incs)), 'Incidence',
                          labels=['%d' % (x + 1) for x in range(max_hh_size + 1)],
                          x_offset=x_offset, legend_title='Household size')

        #output_hh_size_incidence(hh_inc_dat, self.entry_states, 
        #        os.path.join(p['prefix'], 'cases_incidence_by_hh_size.png'))

        ## infant age incidence
        #print "Generating incidence by age (infants) plot"
        #age_dist_quarters = self.get_mean_age_dist_month(times['st'], times['et'])
        #age_inc_dat_quarters = self.get_incidence_by_age_infant(age_dist_quarters, times['st'], times['et'], 2)

        #incs_snapshots = [age_inc_dat_quarters[self.entry_state]]

        #incs_snapshots = self.get_incidence_by_age_infant_snapshots(times['st'], times['et'], 10)

        #plot_statistics(os.path.join(p['prefix'], 'case_incidence_by_age_infant_snapshots.png'),
        #                range(len(incs_snapshots[0])), incs_snapshots, 'Age (months)', 'Incidence',
        #                series_labels=['t=%d' % (x_offset + (x / p['t_per_year'])) for x in t_snapshots],
        #                legend_loc='upper left')

        ## clustering of household susceptibility
        #print "Generating plot for clustering of household susceptibility"

#        hh_at_risk = self.get_series(self.get_hh_at_risk_clustering, times['st'], times['et'], times['t_per_year'] * 5)
#        hh_at_risk_prop = self.get_series(self.get_hh_at_risk_prop, times['st'], times['et'], times['t_per_year'] * 10)

        pts_per_year = 0.1

        # output household-attributable fractions
        print("Generating hh fraction plot")
        hh_fracs, hh_fracs_err = list(zip(*self.get_series(
            self.get_avg_hh_frac, times['st'], times['et'], int(p['t_per_year'] / pts_per_year))))
        output_timeseries(os.path.join(p['prefix'], 'cases_hh_fracs.png'), times,
                          hh_fracs, 'Household fraction',
                          x_offset=x_offset)  # , y_errs=hh_fracs_err)

        print("Generating hh fraction by hh size plot")
        values = self.get_series_dict(self.get_hh_hh_frac, times['st'], times['et'], int(times['t_per_year'] / pts_per_year))
        series = [values[x] for x in sorted(values)]
        labels = ['%d' % x for x in sorted(values)]
        labels[-1] += '+'
        output_timeseries(os.path.join(p['prefix'], 'cases_hh_frac_hh.png'),
                          times, series, 'Household fraction', labels=labels,
                          x_offset=x_offset)

        # output average age of first infection
        print("Generating avg age at first infection plot")
        avg_age, avg_age_err = list(zip(*self.get_series(
            self.get_avg_age_first_infection, times['st'], times['et'], int(times['t_per_year'] / pts_per_year))))
        output_timeseries(os.path.join(p['prefix'], 'cases_age_first_infection.png'), times,
                          avg_age, 'Age',
                          x_offset=x_offset)  # , y_errs=avg_age_err)

        print("Generating avg age at first infection by hh size plot")
        values = self.get_series_dict(
            self.get_hh_avg_age_first_infection, times['st'], times['et'], int(times['t_per_year'] / pts_per_year))
        series = [values[x] for x in sorted(values)]
        labels = ['%d' % x for x in sorted(values)]
        labels[-1] += '+'
        output_timeseries(os.path.join(p['prefix'], 'cases_age_first_infection_hh.png'),
                          times, series, 'Age', labels=labels,
                          x_offset=x_offset)

        
