"""
Observer object for disease state broken down by age and household type
(at the moment, simply with/without children).

cases observer can do this for incidence, but not for prevalence, immunity or susceptibility.

Needs updating: (NON-FUNCTIONAL)
    * MAJOR: snapshots not currently stored in hd5 format!!!
    * implement output_all function
    
"""
from collections import defaultdict
import numpy as np
from disease.observers.obs_base import Observer
import tables as tb


class Case(tb.IsDescription):
    t = tb.UInt32Col()
    
    
class AgeDiseaseObserver(Observer):
    def __init__(self, h5file, t_per_year, t_step):
        super(AgeDiseaseObserver, self).__init__(h5file, 'age_disease', Snapshot, 'State Snapshots')
        self.t_step = t_step        # frequency (in days) with which to take snapshot
        self.t_per_year = t_per_year
        self.snapshots = []

    def create_storage(self, description, title):
        """
        Called when initialised with a h5file opened in 'w'rite mode.
        Adds tables for storage of data (plus local links)
        """
        
        Observer.create_storage(self, description, title)
        group = self.h5file.get_node('/', 'age_disease')
        self.age_dists = self.h5file.create_table(group, 'age_dists',
                    {'t':tb.UInt32Col(), 'dist':tb.UInt32Col(shape=(101))}, 'Age Distributions')
        self.age_dist = self.age_dists.row
        self.age_dists_month = self.h5file.create_table(group, 'age_dists_month',
                    {'t':tb.UInt32Col(), 'dist':tb.UInt32Col(shape=(20))}, 'Age Distributions (months)')
        self.age_dist_month = self.age_dists_month.row
        self.hh_size_dists = self.h5file.create_table(group, 'hh_size_dists',
                    {'t':tb.UInt32Col(), 'dist':tb.UInt32Col(shape=(20))}, 'Household Size Distributions')
        self.hh_size_dist = self.hh_size_dists.row        
        
    def load_storage(self):
        """
        Called when initialised with a h5file opened in 'r'ead or 'a'ppend mode.
        Creates local links to existing tables.
        """
       
        Observer.load_storage(self) 
        self.age_dists = self.h5file.root.cases.age_dists
        self.age_dist = self.age_dists.row        
        self.age_dists_month = self.h5file.root.cases.age_dists_month
        self.age_dist_month = self.age_dists_month.row
        self.hh_size_dists = self.h5file.root.cases.hh_size_dists
        self.hh_size_dist = self.hh_size_dists.row
        
        
    def update(self, t, disease, pop, **kwargs):
        if (t%self.t_step != 0):
            return

        snapshot = []
        for ind in list(pop.I.values()):
            snapshot.append({
                'ID': ind.ID,
                'state': ind.state,
                'age': ind.age,
                'hh_size': pop.hh_size(ind),
                'hh_type': pop.get_hh_type(ind)
                })

        self.snapshots.append((t, snapshot))
        
        self.age_dist['t'] = t 
        self.age_dist['dist'] = pop.age_dist(norm=False)[0]
        self.age_dist.append()
        self.age_dist_month['t'] = t
        self.age_dist_month['dist'] = pop.age_dist_month(bin_days=91, max_age=5, norm=False)[0]
        self.age_dist_month.append()
        self.hh_size_dist['t'] = t 
        self.hh_size_dist['dist'] = pop.group_size_dist('household', 20, norm=False)[0]
        self.hh_size_dist.append()
            
        self.h5file.flush()



    def get_age_dist_bin(self, cur_t, age_bins):
        raw_dist = self.age_dists.where('t==cur_t')
        bin_dist = []
        for cur_bin_l, cur_bin_r in zip(age_bins, age_bins[1:]+[101]):
            bin_dist.append(sum(raw_dist['dist'][cur_bin_l:cur_bin_r]))
        return bin_dist


    def get_state_props_by_age_time(self, state_labels, age_bins):

        props = []
        for ss in self.snapshots:
            cur_t = ss[0]
            cur_dist = self.get_age_dist_bin(cur_t, age_bins)
            d, b = np.histogram(
                [x['age'] for x in ss[1] if x['state'].label in state_labels],
                bins=age_bins, normed=False)
            props.append([float(x)/y for x, y in zip(d, cur_dist)])
        return props


    def get_state_props_by_hh_size_time(self, state_labels):
        
        props = []
        for ss in self.snapshots:
            cur_t = ss[0]
            cur_dist = self.hh_size_dists.where('t==cur_t')
            d, b = np.histogram(
                [x['hh_size'] for x in ss[1] if x['state'].label in state_labels],
                normed = False)
            props.append([float(x)/y for x, y in zip(d, cur_dist)])
        return props
        
    
    def get_avg_age_state_by_hh_size_time(self, state_labels):
        """
        Get the average age of an individual in the specified state(s) by household size over time
        """
        
        ages_raw = []
        for ss in self.snapshots:
            cur_ages = defaultdict(list)
            for cur_ind in list(ss[1].values()):
                if cur_ind['state'] in state_labels:
                    cur_ages[cur_ind['hh_size']].append(cur_ind['age'])
            ages_raw.append(cur_ages)
            
        ages_mean = []
        ages_std = []
        for cur_size in range(1, 7):
            ages_mean.append([np.mean(x.get(cur_size, 0)) for x in ages_raw])
            ages_std.append([np.std(x.get(cur_size, 0)) for x in ages_raw])
        
        return ages_mean, ages_std
        
    


    ### old code (from disease_experiments.output_single) for writing data from this observer.
    """
    # age_disease hasn't been updated to implement output_all yet, as
    # still need to resolve dependency on cases for age_distribution...
    if 'age_disease' in disease.observers:
        # following lines collapse an age distribution into coarser bins.
        # If there isn't a pre-existing function for this, should be moved
        # into a processing library
        print "Processing age_disease output"
        age_dist = []
        new_age_classes = [0, 5, 10, 20, 40, 101]
        cur_sum = 0
        cur_age_class = 1

        for cur_age, count in enumerate(disease.observers['cases'].get_mean_age_dist(st, et)):
            if cur_age >= new_age_classes[cur_age_class]:
                cur_age_class += 1
                age_dist.append(cur_sum)
                cur_sum = 0
            cur_sum += count
        age_dist.append(cur_sum)

        output_age_inc_over_time(
                disease.observers['age_disease'].get_state_series_by_age(
                    ('I'),#,'Ie'),
                    new_age_classes,
                    age_dist),
                new_age_classes,
                times['si'], times['ei'], times['st'], times['et'], times['years_per_tick'], 'Incidence (per 100,000)',
                os.path.join(p['prefix'], 'age_prevalence.png'))#, ymax=0.0002)

        output_age_inc_over_time(
                disease.observers['age_disease'].get_state_series_by_age(
                    ('S','Se'),
                    new_age_classes,
                    age_dist),
                new_age_classes,
                times['si'], times['ei'], times['st'], times['et'], times['years_per_tick'], 'Susceptibility (per 100,000)',
                os.path.join(p['prefix'], 'age_susceptibility.png'))#, ymax=0.009)
        """
        
#      output_state_snapshot(disease.observers['age_disease'].snapshots[0][1],
#              disease.state_labels(), range(0,101,5),
#              os.path.join(p['prefix'], 'age_disease_snapshot.png'))