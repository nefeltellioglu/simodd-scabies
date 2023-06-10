"""
Observer for status (immune, susceptible, etc.) by household size
and other condition (e.g., age, sex, etc.)
"""
import tables as tb
import os

from disease.observers.obs_base import Observer
from disease.experiments.output_disease import *


def always_true(x):
    return True


def not_at_risk(x):
    return not x.state.at_risk


def prev_inf(x):
    return len(x.infections) > 0


def prev_vac(x):
    return len(x.vaccines_received) > 0


def poss_mother(x):
    return (18 < x.age) < 40 and (x.sex == 0)


class HHStatusObserver(Observer):
    """
    """

    def __init__(self, h5file, label='hh_status', base_condition=None,
                 status_condition=None, interval=7, max_hh_size=7):
        self.interval = interval  # recording interval (passed in by days)
        self.max_hh_size = max_hh_size  # max household size (larger households are aggregated)
        self.base_condition = base_condition  # counted in denominator
        self.status_condition = status_condition  # counted in numerator
        desc = {
            't': tb.UInt32Col(),
            'num': tb.UInt32Col(shape=(self.max_hh_size,)),
            'denom': tb.UInt32Col(shape=(self.max_hh_size,))}
        super(HHStatusObserver, self).__init__(h5file, label, desc,
                                               'HH Status Observer: %s' % label)

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
        by_size = pop.all_individuals_by_group_size('household', 7)
        # num_denom is a tuple of (numerator, denominator) tuples, ordered by household size
        num_denom = list((self.calc_matches(by_size[x]) for x in sorted(by_size)))
        while len(num_denom) < 7:
            num_denom.append((0, 0))
        self.row['t'] = t
        self.row['num'] = list(zip(*num_denom))[0]
        self.row['denom'] = list(zip(*num_denom))[1]
        self.row.append()
        self.h5file.flush()

    def get_pos_props(self, window=None):
        props = [x['num'].astype(float) / x['denom'] for x in self.data]
        props = list(zip(*props))
        if window is not None:
            avg_props = []
            for cur_props in props:
                num_points = len(cur_props) / window
                avg_props.append([np.mean(cur_props[x * window:y * window])
                                  for x, y in zip(list(range(0, num_points)), list(range(1, num_points + 1)))])
            return avg_props
        return props


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

        x_offset = 1910-(times['st'] / p['t_per_year'])

        output_timeseries(os.path.join(p['prefix'], '%s_prop.png' % self.label),
                          times, props, 'Proportion immune',
                          labels=['%d' % x for x in range(1, self.max_hh_size)] + ['%d+' % self.max_hh_size],
                          x_offset=x_offset)


        denoms = list(zip(*self.get_pos_denom()))
        output_timeseries(os.path.join(p['prefix'], '%s_prop_denom.png' % self.label),
                          times, denoms, 'Number of people',
                          labels=['%d' % x for x in range(1, self.max_hh_size)] + ['%d+' % self.max_hh_size],
                          x_offset=x_offset)

        #hh_size_props = zip(*self.get_denom_props())

        #output_stacked_bars_timeseries(
        #    os.path.join(p['prefix'], 'demog_hh_size_people_props.png'), times,
        #    hh_size_props, 'Proportion of population',
        #    labels=['%d' % x for x in range(1, self.max_hh_size)] + ['%d+' % self.max_hh_size])
