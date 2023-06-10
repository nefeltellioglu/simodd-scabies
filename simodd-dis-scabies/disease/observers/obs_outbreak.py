import tables as tb

"""
Observer for outbreaks.
"""
from disease.observers.obs_base import Observer

class OutbreakObserver(Observer):
    def __init__(self, h5file, threshold=0):
        self.threshold = threshold
        self.outbreak_present = False       # are we currently recording an outbreak?
        desc = {}
        desc['t_begin'] = tb.UInt32Col()
        desc['t_end'] = tb.UInt32Col()
        desc['t_peak'] = tb.UInt32Col()
        desc['peak_case_count'] = tb.UInt32Col()
        desc['total_case_count'] = tb.UInt32Col()
        super(OutbreakObserver, self).__init__(h5file, 'outbreak', desc, 'Outbreak Observer')
        

    def update(self, t, pop, disease, cases, **kwargs):
        # NB: using cases, rather than new_I, because want to capture outbreaks beginning with
        # BUT: may need to look more closely at this if using an SEIR model
        # first transmission after a case importation
        
        # if disease currently present in population (at sufficient threshold)
        if disease.states['I'].count > self.threshold:
            # and we are not currently recording an outbreak
            if not self.outbreak_present:
                # begin of new outbreak
                self.outbreak_present = True
                self.row['t_begin'] = t
                self.row['peak_case_count'] = 0
                self.row['total_case_count'] = 0
                
            # add new cases (NB: incident, not prevalent!)
            self.row['total_case_count'] += len(cases)
            
            # is this a peak of prevalence?
            if disease.states['I'].count > self.row['peak_case_count']:
                self.row['peak_case_count'] = disease.states['I'].count
                self.row['t_peak'] = t
            

        elif self.outbreak_present:
            # end of current outbreak
            self.outbreak_present = False
            self.row['t_end'] = t
            # don't add 'non-starter' outbreaks
            if self.row['total_case_count'] > 0:
                self.row.append()
                self.h5file.flush()


    def get_outbreak_durations(self):
        """
        Get list of outbreak durations.
        Note that may need to censor beginning to avoid partial outbreaks.
        """
        return self.data.col('t_end') - self.data.col('t_begin')
    

    def get_interval_durations(self):
        """
        Get durations of inter-outbreak intervals.
        """
        ends = self.data.col('t_end')
        begins = self.data.col('t_begin')
        return begins[:-1] - ends[1:]
    
    
    def get_peak_interval_durations(self):
        """
        Get durations of inter-peak intervals.
        """
        peaks = self.data.col('t_peak')
        return peaks[1:] - peaks[:-1]
    
                               
    