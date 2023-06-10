"""
Observer object for disease state (S, E, I, R, etc.) counts.
"""
from disease.observers.obs_base import Observer
import tables as tb
import numpy as np
from . import peakdetect as pd
import os
from disease.experiments.output_disease import output_timeseries
import matplotlib.pyplot as plt


class DiseaseObserver(Observer):
    def __init__(self, h5file, state_labels, disease_states):
        self.state_labels = ['t'] + state_labels
        self.disease_states = disease_states
        desc = dict((x, tb.UInt32Col(pos=i)) for i, x in enumerate(self.state_labels))
        super(DiseaseObserver, self).__init__(h5file, 'disease', desc, 'Disease Observer')

    def update(self, t, disease, **kwargs):
        self.row['t'] = t
        for x in self.state_labels[1:]:
            self.row[x] = disease.states[x].count
        self.row.append()
        self.h5file.flush()
#        if introduction:
#            self.introductions.append((introduction, t))

    def get_counts_by_state(self, label):
        return self.data.col(label)

    def get_props_by_state(self):
        """
        Return a dictionary mapping state labels to proportions.
        """
        sizes = np.sum([self.data.col(x) for x in self.state_labels[1:]],axis=0)
        return dict((x, self.data.col(x)/sizes.astype(float)) for x in self.state_labels[1:])

    def get_pop_sizes(self):
        return np.sum([self.data.col(x) for x in self.state_labels[1:]],axis=0)

    def output_all(self, p, times):
        props = self.get_props_by_state()

        x_offset = 1910-(times['st'] // p['t_per_year'])

        # plot prevalence
        output_timeseries(os.path.join(p['prefix'], 'disease_prevalence.png'),
                          times, [props[x] for x in self.disease_states], 'Prevalence', labels=self.disease_states,
                          x_offset=x_offset)
    
        # plot all disease states as proportions
        output_timeseries(os.path.join(p['prefix'], 'disease_all_states.png'),
                          times, [props[x] for x in self.state_labels[1:]],
                          'Proportion', labels=self.state_labels[1:], logscale=True,
                          x_offset=x_offset)

        if 'Ie' in props:
            props['I'] = np.array(props['In']) + np.array(props['Ie'])

        #print disease.observers['disease'].get_peaks(15, times['st'], times['st'] + times['t_per_year'] * 10)
        #disease.observers['disease'].sweep_peaks(range(10, 160, 10), range(0, 75, 5), times['st'], times['st'] + times['t_per_year'] * 100)

    def get_first_fadeout_time(self, state_labels):
        timeseries = sum(self.get_counts_by_state(y) for y in state_labels)
        try:
            zeros = np.where(timeseries==0)[0]
            #print len(zeros)
            if len(zeros) > 0:
                first_zero = zeros[0]
            else:
                first_zero = len(timeseries)-1
        except Exception:
            print("first zero problem!")
        time = self.data[first_zero]['t']
        return time
    
    def get_max(self, state_label):
        timeseries = [x for x in self.get_counts_by_state(state_label)]
        return max(timeseries)

    def get_time_at_max(self, state_label):
        timeseries = [x for x in self.get_counts_by_state(state_label)]
        peak_inc = max(timeseries)
        index = timeseries.index(peak_inc)
        time = self.data[index]['t']
        return time

    ### Experimental peak detection stuff below here...

    def get_peaks(self, lookahead, delta, st, et):
        # get times and prevalence counts
        x_val = [x['t'] for x in self.data.where('st <= t < et')]
        y_val = [x['I'] for x in self.data.where('st <= t < et')]
        
        _max, _min = pd.peakdetect(y_val, x_val, lookahead, delta)
        
        peak_times, peak_counts = list(zip(*_max))
        peak_intervals = [x-y for x,y in zip(peak_times[1:], peak_times[:-1])]
        
        num_peaks = len(peak_counts)
        mean_interval, stdev_interval = np.mean(peak_intervals), np.std(peak_intervals)
        mean_count, stdev_count = np.mean(peak_counts), np.std(peak_counts)
        
        return num_peaks, (mean_interval, stdev_interval), (mean_count, stdev_count)
    
        
    def sweep_peaks(self, la_vals, delta_vals, st, et):
        
        num_peaks_a = []
        mean_interval_a = []
        std_interval_a = []
        mean_count_a = []
        std_count_a = []
        
        for cur_la in la_vals:
            cur_peaks = []
            cur_interval = []
            cur_interval_std = []
            cur_count = []
            cur_count_std = []
            for cur_delta in delta_vals:
                cur_results = self.get_peaks(cur_la, cur_delta, st, et)
                cur_peaks.append(cur_results[0])
                cur_interval.append(cur_results[1][0])
                cur_interval_std.append(cur_results[1][1])
                cur_count.append(cur_results[2][0])
                cur_count_std.append(cur_results[2][1])
            num_peaks_a.append(cur_peaks)
            mean_interval_a.append(cur_interval)
            std_interval_a.append(cur_interval_std)
            mean_count_a.append(cur_count)
            std_count_a.append(cur_count_std)
        
        fig=plt.figure()
        
        ax = fig.add_subplot(111)
        img = ax.pcolor(np.array(num_peaks_a), cmap='jet')
        ax.set_title('Number of outbreaks')
        ax.set_xlabel('delta')
        ax.set_ylabel('lookahead')
        ax.set_xticks(list(range(len(delta_vals))))
        ax.set_yticks(list(range(len(la_vals))))
        ax.set_xticklabels(delta_vals)
        ax.set_yticklabels(la_vals)
        fig.colorbar(img)
        fig.savefig('num_peaks.png')
        
        fig=plt.figure()
        ax = fig.add_subplot(111)
        img = ax.pcolor(np.array(mean_interval_a), cmap='jet')
        ax.set_title('Mean interval (weeks)')
        ax.set_xlabel('delta')
        ax.set_ylabel('lookahead')
        ax.set_xticks(list(range(len(delta_vals))))
        ax.set_yticks(list(range(len(la_vals))))
        ax.set_xticklabels(delta_vals)
        ax.set_yticklabels(la_vals)
        fig.colorbar(img)
        fig.savefig('interval_mean.png')
        
        fig=plt.figure()
        ax = fig.add_subplot(111)
        img = ax.pcolor(np.array(mean_count_a), cmap='jet')
        ax.set_title('Mean peak prevalence (cases)')
        ax.set_xlabel('delta')
        ax.set_ylabel('lookahead')
        ax.set_xticks(list(range(len(delta_vals))))
        ax.set_yticks(list(range(len(la_vals))))
        ax.set_xticklabels(delta_vals)
        ax.set_yticklabels(la_vals)
        fig.colorbar(img)
        fig.savefig('count_mean.png')
        
        fig=plt.figure()
        ax = fig.add_subplot(111)
        img = ax.pcolor(np.array(std_interval_a), cmap='jet')
        ax.set_title('SD interval (weeks)')
        ax.set_xlabel('delta')
        ax.set_ylabel('lookahead')
        ax.set_xticks(list(range(len(delta_vals))))
        ax.set_yticks(list(range(len(la_vals))))
        ax.set_xticklabels(delta_vals)
        ax.set_yticklabels(la_vals)
        fig.colorbar(img)
        fig.savefig('interval_std.png')
        
        fig=plt.figure()
        ax = fig.add_subplot(111)
        img = ax.pcolor(np.array(std_count_a), cmap='jet')
        ax.set_title('SD peak prevalence (cases)')
        ax.set_xlabel('delta')
        ax.set_ylabel('lookahead')
        ax.set_xticks(list(range(len(delta_vals))))
        ax.set_yticks(list(range(len(la_vals))))
        ax.set_xticklabels(delta_vals)
        ax.set_yticklabels(la_vals)
        fig.colorbar(img)
        fig.savefig('count_std.png')            
        
    
        



#    def get_peaks(self, threshold, st, et):
#        """
#        Get the times of each peak prevalence between which prevalence drops below the specified cutoff.
#        """
#        peaks = []
#        times = []
#        
#        outbreak_present = False
#        
#        #peak_count = 0
#        #peak_time = 0
#        
#        for cur in self.data.where('st <= t < et'):
#            if cur['I'] > threshold:
#                if not outbreak_present:
#                    outbreak_present = True
#                    peak_count = 0
#                
#                if cur['I'] > peak_count:
#                    peak_count = cur['I']
#                    peak_time = cur['t']
#                    
#            elif outbreak_present:
#                outbreak_present = False
#                peaks.append(peak_count)
#                times.append(peak_time)
#                
#        return peaks, times
            
    