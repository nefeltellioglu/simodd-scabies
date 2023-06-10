import tables as tb
import os

import numpy as np

from disease.observers.obs_base import Observer
from disease.experiments.output_disease import *


class HHSusceptibleObserver(Observer):
    """
    Observer for tracking number of susceptible individuals by household size
    """

    def __init__(self, h5file, interval=1, max_size=7):
        self.interval = interval  # recording interval
        self.max_size = max_size  # maximum household size to track
        desc = {'t': tb.UInt32Col()}  # time of observation
        # number of households of size x (for x in range 1...max_size)
        desc.update(dict(('count_%d' % x, tb.UInt32Col()) for x in range(1, max_size + 1)))
        # vector of susceptible counts for households of size x
        desc.update(dict(('susc_%d' % x, tb.UInt32Col(shape=(x + 1))) for x in range(1, max_size + 1)))
        desc.update(dict(('age_mean_%d' % x, tb.Float32Col()) for x in range(1, max_size + 1)))
        desc.update(dict(('age_std_%d' % x, tb.Float32Col()) for x in range(1, max_size + 1)))
        super(HHSusceptibleObserver, self).__init__(h5file, 'hh_suscept', desc, 'Household susceptibility')


    def update(self, t, pop, **kwargs):
        if t % self.interval > 0: return

        self.row['t'] = t
        hh_by_size = pop.all_groups_by_size('household', self.max_size)
        for cur_size, cur_hhs in list(hh_by_size.items()):
            susc_ages = []
            self.row['count_%d' % cur_size] = len(cur_hhs)
            susc_count = np.zeros((cur_size + 1), dtype=np.int)
            for cur_hh in cur_hhs:
                susc_ages.extend([x.age for x in cur_hh if x.state.at_risk])
                cur_count = min(len([x for x in cur_hh if x.state.at_risk]), self.max_size)
                #print cur_size, len(cur_hh), cur_count
                susc_count[cur_count] += 1
            #print susc_count, type(susc_count)
            self.row['susc_%d' % cur_size] = susc_count.tolist()
            self.row['age_mean_%d' % cur_size] = np.mean(susc_ages)
            self.row['age_std_%d' % cur_size] = np.std(susc_ages)
        self.row.append()
        self.h5file.flush()


    def get_susceptible_props(self, hh_size, avg=None):
        # returns a vector of the proportions of households of size 'hh_size' with n susceptibles
        # avg is number of data points to average over
        counts = {}
        for num_susc in range(0, hh_size + 1):
            counts[num_susc] = [float(x[num_susc]) / (y) for x, y in
                                zip(self.data.col('susc_%d' % hh_size), self.data.col('count_%d' % hh_size))]
        if avg is not None:
            avg_counts = {}
            for k, v in counts:
                num_points = len(v) / avg
                avg_counts[k] = [np.mean(v[x * avg:y * avg]) for x, y in
                                 zip(list(range(0, num_points)), list(range(1, num_points + 1)))]
            return avg_counts
        return counts


    def get_total_susceptible_props(self):
        # returns the total proportion of susceptible individuals in households of size 'hh_size'
        susc_counts = {}
        total_counts = {}
        susc_props = {}
        for cur_hh_size in range(1, self.max_size + 1):
            susc_counts[cur_hh_size] = np.zeros(len(self.data.col('t')))
            total_counts[cur_hh_size] = np.array([cur_hh_size * x for x in self.data.col('count_%d' % cur_hh_size)])
            for num_susc in range(0, cur_hh_size + 1):
                susc_counts[cur_hh_size] += [num_susc * x[num_susc] for x in self.data.col('susc_%d' % cur_hh_size)]
            susc_props[cur_hh_size] = (susc_counts[cur_hh_size] / total_counts[cur_hh_size]).tolist()
        return susc_props


    def output_all(self, p, times):
        """
        Write a series of plots showing susceptibility proportions over time, by household size.
        """

        for hh_size in range(1, self.max_size + 1):
            props = self.get_susceptible_props(hh_size)
            labels = ['%d' % x for x in sorted(props.keys())]
            values = [props[x] for x in sorted(props.keys())]
            output_timeseries(os.path.join(p['prefix'], 'hh_susc_props_%d.png' % hh_size), times,
                              values, 'Proportion', labels=labels, logscale=True,
                              title='Proportion of %d-person households with n susceptibles' % hh_size)

        total_props = self.get_total_susceptible_props()
        labels = ['%d' % x for x in sorted(total_props.keys())]
        values = [total_props[x] for x in sorted(total_props.keys())]
        output_timeseries(os.path.join(p['prefix'], 'hh_total_susc_props.png'), times,
                          values, 'Proportion susceptible', labels=labels, logscale=True,
                          title='Proportion of people susceptible by household size')

        if 'age_mean_1' in self.data:
            output_timeseries(os.path.join(p['prefix'], 'hh_susc_ages.png'), times,
                              [self.data.col('age_mean_%d' % hh_size).tolist() for hh_size in
                               range(1, self.max_size + 1)],
                              'Age',
                              #y_errs=[self.data.col('age_std_%d'%hh_size).tolist() for hh_size in range(1, self.max_size+1)],
                              labels=['%d' % x for x in range(1, self.max_size + 1)],
                              title='Mean age of susceptible individuals by household size',
            )
            