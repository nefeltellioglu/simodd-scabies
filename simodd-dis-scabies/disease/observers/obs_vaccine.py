"""
Observer for vaccine state snapshots.
"""
import os
from disease.observers.obs_base import Observer
from disease.experiments.output_disease import output_timeseries
import tables as tb
import numpy as np

class VaccineObserver(Observer):
    def __init__(self, h5file, vac, interval):
        self.vac_labels = [x.label for x in vac]
        self.interval = interval
        desc = dict((x, tb.UInt32Col()) for x in self.vac_labels)
        desc.update(dict(('%s_new' % x, tb.UInt32Col()) for x in self.vac_labels))
        desc.update(dict(('%s_ever' % x, tb.UInt32Col()) for x in self.vac_labels))
        desc.update(dict(('%s_eff' % x, tb.UInt32Col()) for x in self.vac_labels))
        desc['t'] = tb.UInt32Col()
        desc['birth_count'] = tb.UInt32Col()
        desc['pop_size'] = tb.UInt32Col()
        super(VaccineObserver, self).__init__(h5file, 'vaccine', desc, 'Vaccine Observer')
    
        
    def update(self, t, pop, disease, **kwargs):
        if t % self.interval > 0: return
        self.row['t'] = t
        for v in list(disease.vaccines.values()):
            #recpt = [x.ID for x in pop.I.values() if v.label in x.vaccines_received]
            cov = len([x for x in list(pop.I.values()) if v.label in x.vaccines_received])
            self.row[v.label] = cov
            eff = len([x for x in list(pop.I.values()) if v.label in x.vaccines_received
                       and x.vaccines_effective[x.vaccines_received.index(v.label)]])
            self.row['%s_eff'%v.label] = eff
            self.row['%s_new'%v.label] = v.count
            self.row['%s_ever'%v.label] = v.count_ever
        self.row['pop_size'] = len(pop.I)
        self.row['birth_count'] = disease.birth_count   # NB: this only captures instantaneous birth count...
        self.row.append()
        self.h5file.flush()
        
    
    def get_cov_props(self, label):
        return self.data.col(label).astype(np.float) / self.data.col('pop_size')
    
    
    def get_birth_props(self, label):
        return [(float(x)/y if y > 0 else np.nan) \
            for x, y in zip(self.data.col('%s_new'%label), self.data.col('birth_count'))]
        

    def output_all(self, p, times):
    
        vacc_cov = []
        new_vacc = []
        for cur_v in self.vac_labels:
            vacc_cov.append(self.get_cov_props(cur_v))
            
            #new_vacc.append([cur_v, self.get_birth_props(cur_v)])
            #new_vacc = []

        x_offset = 1710

        output_timeseries(os.path.join(p['prefix'], 'vaccine_coverage.png'), times, vacc_cov,
                          'Coverage', labels=self.vac_labels, x_offset=x_offset)
