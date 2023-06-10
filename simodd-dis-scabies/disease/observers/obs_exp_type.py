"""
Observer for information on exposure types.
"""
import numpy as np
from disease.observers.obs_base import Observer
from collections import defaultdict
import itertools
from pylab import NaN
import time


class ExpTypeObserver(Observer):
    def __init__(self):
        super(ExpTypeObserver, self).__init__('exp_type')
        self.data = {'infection_kids':0, 'infection_nokids':0,
                     'boosting_kids':0, 'boosting_nokids':0}


    def process_case(self, t, pop, ind, exp_type):
        if ind.age > 17:
            hh_type = 'kids' if pop.get_hh_type(ind).endswith('kids') else 'nokids'
            exp_key = exp_type + '_' + hh_type
#        print exp_key
            self.data[exp_key] += 1


    def update(self, t, pop, cases, boosting, **kwargs):
#        print "%%%", cases
        for (ind, hh_frac) in cases:
            self.process_case(t, pop, ind, 'infection')
        for (ind, hh_frac) in boosting:
            self.process_case(t, pop, ind, 'boosting')


    def reset(self):
        """Remove all stats"""
        self.data = []

    def output_all(self, p, times):
        self.write_ratios(os.path.join(p['prefix'], 'exp_type_ratios.csv'))
    

    def write_ratios(self, filename):
        ofile = open(filename, 'w')
        ofile.write('%d,%d\n' % (self.data['infection_kids'], self.data['infection_nokids']))
        ofile.write('%d,%d\n' % (self.data['boosting_kids'], self.data['boosting_nokids']))
        ofile.write('%g,%g\n' % (float(self.data['infection_kids']) / (self.data['infection_kids'] + self.data['boosting_kids']),
                                 float(self.data['infection_nokids']) / (self.data['infection_nokids'] + self.data['boosting_nokids'])))
        ofile.close()

 
