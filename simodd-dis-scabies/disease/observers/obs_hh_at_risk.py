from disease.observers.obs_base import Observer
from disease.experiments.output_disease import *
import tables as tb
import numpy as np
import os
from collections import defaultdict

class HHAtRiskObserver(Observer):
    """
    """

    def __init__(self, h5file, interval=7):
        self.interval = interval                    # recording interval
        desc = {
            't': tb.UInt32Col(),
            'hh_id': tb.UInt32Col(),
            'hh_age': tb.UInt32Col(),
            'count': tb.UInt8Col()}   # number of hh members at risk
        super(HHAtRiskObserver, self).__init__(h5file, 'hh_at_risk', desc, 'Household members at risk')


    def update(self, t, pop, **kwargs):
        if t % self.interval > 0: return
        
        for hh_id, cur_hh in list(pop.groups['household'].items()):
            self.row['t'] = t
            self.row['hh_id'] = hh_id
            self.row['hh_age'] = t - pop.households[hh_id].founded
            self.row['count'] = len([x for x in cur_hh if x.state.at_risk])
            self.row.append()
        self.h5file.flush()
            
        
    def get_counts_by_age(self):
        c_by_age = defaultdict(list)
        for x in self.data:
#            if x['count'] > 0:
            c_by_age[x['hh_age'] / 364].append(x['count'])
        return c_by_age
        
        
    def get_counts_by_hh(self):
        c_by_hh = defaultdict(list)
        for x in self.data:
            c_by_hh[x['hh_id']].append(x['count'])
        return c_by_hh
        
        
    def output_all(self, p, times):
        output_hh_risk(self.get_counts_by_age(), os.path.join(p['prefix'], 'hh_risk.png'))
        output_hh_risk_series(self.get_counts_by_hh(), os.path.join(p['prefix'], 'hh_risk_series.png'))       
        
        
        
#        by_size = pop.all_individuals_by_group_size('household', 7)
#        # num_denom is a tuple of (num, denom) tuples, ordered by household size
#        num_denom = list((self.calc_matches(by_size[x]) for x in sorted(by_size)))
#        self.snapshot['t'] = t
#        self.snapshot['num'] = zip(*num_denom)[0]
#        self.snapshot['denom'] = zip(*num_denom)[1]
#        self.snapshot.append()
#        self.h5file.flush()
