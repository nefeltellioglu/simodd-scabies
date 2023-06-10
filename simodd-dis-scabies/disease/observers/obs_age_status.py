import tables as tb

import numpy as np
import os

from disease.observers.obs_base import Observer
from disease.experiments.output_disease import *


def always_true(x):
    return True


def not_at_risk(x):
    return not x.state.at_risk


class AgeStatusObserver(Observer):
    """
    """

    def __init__(self, h5file, label='hh_status', base_condition=None,
                 status_condition=None, interval=7, age_bins=None):
        self.interval = interval  # recording interval (passed in by days)
        self.age_bins = age_bins  # max household size (larger households are aggregated)
        self.base_condition = base_condition  # counted in denominator
        self.status_condition = status_condition  # counted in numerator
        desc = {
            't': tb.UInt32Col(),
            'num': tb.UInt32Col(shape=(len(age_bins),)),
            'denom': tb.UInt32Col(shape=(len(age_bins),))}
        super(AgeStatusObserver, self).__init__(
            h5file, label, desc, 'Age Status Observer: %s' % label)

    def calc_matches(self, raw_inds):
        """
        Return the numerator (number of individuals who satisfy both the base condition
        AND the status condition) and denominator (total number of individuals who satisfy
        the base condition).
        """
        inds = [x for x in raw_inds if self.base_condition(x)]
        return len([x for x in inds if self.status_condition(x)]), len(inds)

    def update(self, t, pop, **kwargs):
        if t % self.interval > 0:
            return
        by_size = dict([(x, pop.individuals_by_age(x, y)) for x, y in zip(self.age_bins, self.age_bins[1:] + [102])])
        # num_denom is a tuple of (numerator, denominator) tuples, ordered by household size
        num_denom = list((self.calc_matches(by_size[x]) for x in sorted(by_size)))
        #while len(num_denom) < 7:
        #    num_denom.append((0, 0))
        self.row['t'] = t
        self.row['num'] = list(zip(*num_denom))[0]
        self.row['denom'] = list(zip(*num_denom))[1]
        self.row.append()
        self.h5file.flush()

    def get_pos_props(self):
        return [x['num'].astype(float) / x['denom'] for x in self.data]

    def get_pos_denom(self):
        return self.data.col('denom')

    def get_denom_props(self):
        props = []
        for cur in self.data.col('denom'):
            total = float(np.sum(cur))
            props.append([x / total for x in cur])
        return props

    def output_all(self, p, times):
        props = list(zip(*self.get_pos_props()))

        age_labels = ['%d-%d' % (x, y) for x, y in zip(self.age_bins, self.age_bins[1:] + [100])]

        output_timeseries(os.path.join(p['prefix'], '%s_prop_age.png' % self.label),
                          times, props, 'Proportion immune', labels=age_labels)

        denoms = list(zip(*self.get_pos_denom()))
        output_timeseries(os.path.join(p['prefix'], '%s_prop_age_denom.png' % self.label),
                          times, denoms, 'Number of people', labels=age_labels)

        #hh_size_props = zip(*self.get_denom_props())
        #
        #output_stacked_bars_timeseries(
        #    os.path.join(p['prefix'], 'demog_age_props.png'), times,
        #    hh_size_props, 'Proportion of population', labels = age_labels)