import tables as tb

"""
Observer for importations.
"""
from disease.observers.obs_base import Observer

class ImportsObserver(Observer):
    def __init__(self, h5file):
        desc = {}
        desc['disease_present'] = tb.BoolCol()
        desc['t'] = tb.UInt32Col()
        desc['hh_size'] = tb.UInt32Col()
        desc['S_hh'] = tb.UInt32Col()       # proportion of household currently susceptible
        desc['S_pop'] = tb.UInt32Col()      # proportion of population currently susceptible
        super(ImportsObserver, self).__init__(h5file, 'imports', desc, 'Importations Observer')
       
        
    def update(self, t, disease, pop, introduction, **kwargs):
        if not introduction: return
        self.row['t'] = t
        self.row['disease_present'] = (disease.states['I'].count > 1)
        hh_members = pop.housemates(introduction)
        self.row['hh_size'] = len(hh_members)+1
        # NB: currently assumes that S is the only susceptible state!
        self.row['S_hh'] = len([x for x in hh_members if x.state.label == 'S'])
        self.row['S_pop'] = float(disease.states['S'].count)
        self.row.append()
        self.h5file.flush()
        
        
    