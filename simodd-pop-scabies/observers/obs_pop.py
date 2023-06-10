"""
Observer for demographic (population) information.
"""

import tables as tb

from population.pop_info import age_dist, group_size_dist, group_size_avg
from observers.obs_base import Observer


class PopulationObserver(Observer):
    def __init__(self, h5file, interval=364):
        self.interval = interval
        self.max_hh = 20
        self.age_dists = None
        self.hh_size_dists = None
        self.hh_size_dists_count = None
        desc = {
            't': tb.UInt32Col(),
            'pop_size': tb.UInt32Col(),
            'hh_count': tb.UInt32Col(),
            'age_mean': tb.Float32Col(),
            'hh_size_mean': tb.Float32Col()
        }
        super(PopulationObserver, self).__init__(h5file, 'population', desc, 'Population Observer')

    def create_storage(self, description, title):
        """
        Called when initialised with a h5file opened in 'w'rite mode.
        Adds tables for storage of data (plus local links)
        """

        Observer.create_storage(self, description, title)
        group = self.h5file.get_node('/', 'population')
        self.age_dists = self.h5file.create_table(
            group, 'age_dists',
            {'t': tb.UInt32Col(), 'dist': tb.UInt32Col(shape=(101))}, 'Age Distributions')
        self.hh_size_dists = self.h5file.create_table(
            group, 'hh_size_dists',
            {'t': tb.Float32Col(), 'dist': tb.UInt32Col(shape=(20))}, 'Household Size Distributions')
        self.hh_size_dists_count = self.h5file.create_table(
            group, 'hh_size_dists_count',
            {'t': tb.UInt32Col(), 'dist': tb.UInt32Col(shape=(20))}, 'Household Size Distributions (count)')

    def load_storage(self):
        """
        Called when initialised with a h5file opened in 'r'ead or 'a'ppend mode.
        Creates local links to existing tables.
        """

        Observer.load_storage(self)
        self.age_dists = self.h5file.root.cases.age_dists
        self.hh_size_dists = self.h5file.root.cases.hh_size_dists
        self.hh_size_dists_count = self.h5file.root.cases.hh_size_dists_count

    def update(self, t, pop, **kwargs):
        if t % self.interval > 0:
            return

        self.age_dists.row['t'] = t
        self.age_dists.row['dist'] = age_dist(pop.I, norm=False)[0]
        self.age_dists.row.append()
        self.hh_size_dists.row['t'] = t
        self.hh_size_dists.row['dist'] = group_size_dist(pop.groups, 'household', self.max_hh)[0]
        self.hh_size_dists.row.append()
        self.hh_size_dists_count.row['t'] = t
        self.hh_size_dists_count.row['dist'] = group_size_dist(pop.groups, 'household', self.max_hh, False)[0]
        self.hh_size_dists_count.row.append()

        # self.fam_types.append(self.P.sum_hh_stats_group())

        self.row['t'] = t
        self.row['pop_size'] = len(pop.I)
        self.row['hh_count'] = len(pop.groups['household'])
        self.row['hh_size_mean'] =  group_size_avg(pop.groups, 'household')

        self.row.append()
        self.h5file.flush()
