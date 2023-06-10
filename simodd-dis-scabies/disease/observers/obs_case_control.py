import numpy as np
from disease.observers.obs_base import Observer
from collections import defaultdict
import itertools

"""
Case condition = adult (age > 17)

Age-matched control for each case.

Exposure = in house containing children.
"""

class CaseControlObserver(Observer):
    def __init__(self):
        super(CaseControlObserver, self).__init__('case_control')
        self.case_ids = set()
        self.control_ids = set()
        self.counts = {}
        self.counts['case_E'] = 0
        self.counts['case_U'] = 0
        self.counts['control_E'] = 0
        self.counts['control_U'] = 0
        self.counts['unmatched_cases'] = 0


    def match_case(self, pop, case_ind, rng):
        case_age = case_ind.age
        possible_controls = [x \
                for x in pop.individuals_by_age(case_age, case_age) \
                if x.state.label in ['S', 'Se'] \
                and x.ID not in self.control_ids]
        if not possible_controls:
#            print "Warning: Failed to find an eligible control!"
            return None
        else:
            return rng.choice(possible_controls)



    def update(self, t, pop, cases, **kwargs):
        # filter for eligible cases (ie, above a certain age, etc)
        eligible_cases = [x[0] for x in cases if x[0].age > 17]
        
        # create list of matched controls
        controls = []
        for case_ind in eligible_cases:
            control = self.match_case(pop, case_ind, kwargs['rng'])
            if control:
                controls.append(control)
#                print "CASE_ID (age):", case_ind.ID, "(%d)" % case_ind.age, \
#                        " : CONTROL_ID (age):", controls[-1].ID, "(%d)" % \
#                        controls[-1].age
            else:
                self.counts['unmatched_cases'] += 1

        # count number of cases and controls that are exposed/unexposed
        for case_ind in eligible_cases:
            self.case_ids.add(case_ind.ID)
            if pop.get_hh_type(case_ind).endswith('kids'):
                self.counts['case_E'] += 1
            else:
                self.counts['case_U'] += 1
        for control_ind in controls:
            self.case_ids.add(control_ind.ID)
            if pop.get_hh_type(control_ind).endswith('kids'):
                self.counts['control_E'] += 1
            else:
                self.counts['control_U'] += 1



    def get_relative_risk(self, t_years):
        # number of exposed and unexposed cases
        C_E = float(self.counts['case_E'])
        C_U = float(self.counts['case_U'])

        # size of exposed and unexposed groups
        ### Think this might be incorrect:
        # should be total size of exposed/unexposed groups, not just 
        # controls
        N_E = C_E + self.counts['control_E']
        N_U = C_U + self.counts['control_U']

        # person years at risk for exposed and unexposed groups:
        #   avg number of people at risk during period, 
        #   multiplied by length of period (in years)
        pyar_E = (N_E + self.counts['control_E']) / 2.0
        pyar_U = (N_U + self.counts['control_U']) / 2.0

        relative_risk = (C_E / C_U) / (pyar_E / pyar_U)
        return relative_risk



    def write_numbers(self, ofilename):
        ofile = open(ofilename, 'w')
        ofile.write(',Exposed,Unexposed,Total\n')
        ofile.write('Case,%d,%d,%d\n' % (self.counts['case_E'], self.counts['case_U'],
                                         self.counts['case_E']+self.counts['case_U']))
        ofile.write('Control,%d,%d,%d\n' % (self.counts['control_E'], self.counts['control_U'],
                                            self.counts['control_E']+self.counts['control_U']))
        ofile.write('Total,%d,%d,%d\n' % (self.counts['case_E']+self.counts['control_E'], self.counts['case_U']+self.counts['control_U'],
                                          self.counts['case_E']+self.counts['control_E']+self.counts['case_U']+self.counts['control_U']))
        ofile.write('PYAR,%d,%d\n' % ((self.counts['case_E']+self.counts['control_E']+self.counts['control_E']) / 2.0,
                                      (self.counts['case_U']+self.counts['control_U']+self.counts['control_U']) / 2.0))
        ofile.close()

