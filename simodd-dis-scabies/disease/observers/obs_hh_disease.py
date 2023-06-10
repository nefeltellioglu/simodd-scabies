"""
Observer object for household disease state (S, E, I, R, etc.) counts.
"""
from disease.observers.obs_base import Observer
from collections import defaultdict

class HHDiseaseObserver(Observer):
    def __init__(self, state_labels):
        super(HHDiseaseObserver, self).__init__('disease')
        self.state_labels = ['t'] + state_labels
        self.disease_states = []


    def update(self, t, disease, pop, **kwargs):
        cur_state = dict([(x, []) for x in self.state_labels])
        for hh in pop.groups['household']:
            hh_state = 'R'
            for i_id in hh:
                if i.state.label == 'S':
                    hh_state = 'S'
                if i.state.label == 'I':
                    cur_state = 'I'
                    break
            cur_state[hh_state].append(len(hh))
        self.disease_states.append([t] + [cur_state[k] for k in \
                sorted(cur_state, key=lambda x: disease.state[x].order)])


    def get_counts(self):
        return list(zip(*self.disease_states))


    def get_counts_by_state(self, label):
        return list(zip(*self.disease_states))[self.state_labels.index(label)]
