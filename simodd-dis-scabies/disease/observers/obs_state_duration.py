__author__ = 'ngeard'

"""
Observer object for state durations.
"""
from disease.observers.obs_base import Observer
import tables as tb
import numpy as np
from . import peakdetect as pd
import os
from disease.experiments.output_disease import output_timeseries
import matplotlib.pyplot as plt


class DurationObserver(Observer):
    def __init__(self, h5file, state_labels, vacc_labels=None):
        self.state_labels = state_labels
        self.vacc_labels = vacc_labels
        desc = {}
        super(DurationObserver, self).__init__(h5file, 'duration', desc, 'Duration Observer')

    def create_storage(self, description, title):

        Observer.create_storage(self, description, title)
        group = self.h5file.get_node('/', 'duration')

        self.state_durations = {}
        self.state_d = {}

        for cur_state in self.state_labels:
            self.state_durations[cur_state] = self.h5file.create_table(
                group, 'dur_%s' % cur_state, {'t': tb.UInt32Col(), 'd': tb.UInt32Col()}, 'Durations')
            self.state_d[cur_state] = self.state_durations[cur_state].row

    def update(self, t, pop, **kwargs):
        for ind in list(pop.I.values()):
            for state, duration in list(ind.counters.items()):
                if duration == ind.durations[state] - 1:
                    self.state_d[state]['t'] = t
                    self.state_d[state]['d'] = ind.durations[state]
                    self.state_d[state].append()
        self.h5file.flush()

    def output_all(self, p, times):
        pass