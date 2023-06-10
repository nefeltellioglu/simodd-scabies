"""
Observer object for summary stats:
-Overall prevalence
-Age and hh size dependent prevalence (age04, age59 age1014,..., hh1,hh2,hh3,...)
-Mean and variance of ages of infected individuals
-Mean and variance of hh size of infected individuals
-Prop of hhs with 0 prevalence
-Prop of hh with 1 case

"""
from disease.observers.obs_base import Observer
import tables as tb
import numpy as np
from . import peakdetect as pd
import os
from disease.experiments.output_disease import output_timeseries
import matplotlib.pyplot as plt
import csv
from math import ceil
from collections import defaultdict
from scipy.optimize import curve_fit
import pandas as pd
from pylab import NaN

from population.pop_info import hh_size, age_dist, age_dist_month


class SummaryStatsObserver(Observer):
    """
    An observer for recording summary stats.
    """
    def __init__(self,h5file,observefrom,et,fname,  t_per_year=364):
        self.fname=fname
        self.observefrom=observefrom
        self.et=et
        self.age_dists_prevalence = None
        self.age_dist_prevalence = None
        self.hh_size_dists_prevalence = None
        self.overallprevalence=None
        self.mean_age_infected=None
        self.variance_age_infected=None
        self.mean_hhsize_infected=None
        self.hh_with_zerocase=None #hh percentage with zero case
        self.hh_with_onecase=None #hh percentage with one case
        self.hhtransmission_percentage=None #percentage of I people with hh transmission
        self.t_per_year = t_per_year
        #self.store_source = store_source
        self.stats_df=pd.DataFrame( )
        self.hh_size_dists_prevalence_df=pd.DataFrame()
        self.age_dists_prevalence_df=pd.DataFrame()
        self.summarystats_df=pd.DataFrame()
        self.hh_size_comm_dists_prevalence_df=pd.DataFrame()
        self.age_comm_dists_prevalence_df=pd.DataFrame()
        self.summarystats_df2=pd.DataFrame()
        desc = {
            'age_dists_prevalence': tb.Float32Col(), #Q.typess?
            'hh_size_dists_prevalence': tb.Float32Col(),
            'stats':tb.Float32Col(),
            #'overallprevalence': tb.Float32Col(),
            #'mean_age_infected': tb.Float32Col(),
            #'variance_age_infected': tb.Float32Col(),
            #'mean_hhsize_infected': tb.Float32Col(),
            #'hh_with_zerocase': tb.Float32Col(),
            #'hh_with_onecase': tb.Float32Col(),
            #'hhtransmission_percentage': tb.Float32Col(),
        }   
        super(SummaryStatsObserver, self).__init__(h5file, 'stats', desc, 'Summary Stats') #??? 

    def create_storage(self, description, title): #Q.Shapes?
        """
        Called when initialised with a h5file opened in 'w'rite mode.
        Adds tables for storage of data (plus local links)
        """

        Observer.create_storage(self, description, title)
        group = self.h5file.get_node('/', 'stats')
        
        self.age_dists_prevalence = self.h5file.create_table(group, 'age_dists_prevalence',
                                                 {'t': tb.UInt32Col(), 'dist': tb.Float32Col(shape=(16))},
                                                 'Age Specific Prevalence')
        self.age_dist_prevalence =self.age_dists_prevalence.row
        
        self.hh_size_dists_prevalence = self.h5file.create_table(group, 'hh_size_dists_prevalence',
                                                 {'t': tb.UInt32Col(), 'dist': tb.Float32Col(shape=(12))},
                                                 'HH Size Specific Prevalence')
        self.hh_size_dist_prevalence =self.hh_size_dists_prevalence.row
        
        self.stats = self.h5file.create_table(group, 'stats',
            {'t': tb.UInt32Col(), 
             'overallprevalence': tb.Float32Col(shape=(1)),
            'mean_age_infected': tb.Float32Col(shape=(1)),
             'variance_age_infected':tb.Float32Col(shape=(1)),
             'mean_hhsize_infected':tb.Float32Col(shape=(1)),
             'hh_with_zerocase':tb.Float32Col(shape=(1)),
             'hh_with_onecase':tb.Float32Col(shape=(1)),
             'hh_with_twocase':tb.Float32Col(shape=(1)),
             'com_transmission_percentage':tb.Float32Col(shape=(1)),
            }, 'single Statistics')
        self.stat=self.stats.row
        
        
    def load_storage(self):
        """
        Called when initialised with a h5file opened in 'r'ead or 'a'ppend mode.
        Creates local links to existing tables.
        """
        #Q.row?
        Observer.load_storage(self)
        self.age_dists_prevalence = self.h5file.root.cases.age_dists_prevalence
        self.age_dist_prevalence = self.h5file.root.cases.age_dists_prevalence.row
        self.hh_size_dists_prevalence = self.h5file.root.cases.hh_size_dists_prevalence
        self.hh_size_dist_prevalence = self.h5file.root.cases.hh_size_dists_prevalence.row
        self.stats=self.h5file.root.cases.stats
        self.stat=self.stats.row  
        
    def _reset(self):
        """
        Remove all data and re-create empty storage.
        """
        self.h5file.remove_node('/', 'stats', recursive=True)
        self.create_storage(desc, 'Summary Stats')
        
    def update(self, t, disease, pop,cases, **kwargs):
        """
        Update details of all stats this time step.
        """
        if t % self.t_per_year in range(0, self.t_per_year, self.t_per_year // 4):
            if t >= (self.observefrom):
                #age dependent prevalence
                self.age_dist_prevalence['t'] = t
                ages_infected= np.array([i.age for i in list(pop.I.values()) if i.state.label == "I"])
                ages_infected_hist=np.histogram(ages_infected, bins=16, range=(0, 80))[0]
                ages_all= np.array([i.age for i in list(pop.I.values())])
                ages_all_hist=np.histogram(ages_all, bins=16, range=(0, 80))[0]
                self.age_dist_prevalence['dist'] =np.array([x/y if y else np.nan for x, y in zip(ages_infected_hist, ages_all_hist)])
                
                #display(self.age_dist_prevalence['dist'])
                new_df_ages=pd.DataFrame([[t]]+[[i] for i in self.age_dist_prevalence['dist']]).transpose()
                #display(new_df_ages)
                self.age_dists_prevalence_df=pd.concat([self.age_dists_prevalence_df, new_df_ages], ignore_index=True)
                #display(self.age_dists_prevalence_df)
                self.age_dist_prevalence.append()

                #hh size dependent prevalence
                self.hh_size_dist_prevalence['t'] = t
                hh_size_infected=np.array([hh_size(pop, i) for i in list(pop.I.values()) if i.state.label == "I"])
                hh_size_infected_hist=np.histogram(hh_size_infected, bins=12, range=(1,13))[0]
                hh_size_all= np.array([hh_size(pop, i) for i in list(pop.I.values())])
                hh_size_all_hist=np.histogram(hh_size_all, bins=12, range=(1, 13))[0]
                self.hh_size_dist_prevalence['dist'] =np.array([x/y if y else np.nan for x, y in zip(hh_size_infected_hist, hh_size_all_hist)])
                
                #display(self.hh_size_dist_prevalence['dist'])
                new_df_hhs=pd.DataFrame([[t]]+[[i] for i in self.hh_size_dist_prevalence['dist'] ]).transpose()
                self.hh_size_dists_prevalence_df=pd.concat([self.hh_size_dists_prevalence_df, new_df_hhs], ignore_index=True)
                self.hh_size_dist_prevalence.append()


                #other simple stats
                self.stat['t']=t
                self.stat['overallprevalence']=len([i for i in list(pop.I.values()) if i.state.label == "I"])/len(list(pop.I.values()))
                ages_infected=np.array([i.age for i in list(pop.I.values()) if i.state.label == "I"])
                mean_age_infected=np.mean(ages_infected) if len(ages_infected)>0 else 0
                var_age_infected=np.var(ages_infected) if len(ages_infected)>0 else 0
                self.stat['mean_age_infected']=mean_age_infected
                self.stat['variance_age_infected']=var_age_infected
                hhsizes_infected=np.array([hh_size(pop, i) for i in list(pop.I.values()) if i.state.label == "I"])
                mean_hhsize_infected=np.mean(hhsizes_infected) if len(hhsizes_infected)>0 else 0
                self.stat['mean_hhsize_infected']=mean_hhsize_infected
                
                #comm transmission
                inf_pop=[i for i in list(pop.I.values()) if i.state.label == "I"]
                com_transmisssion_perc=len([i for i in list(pop.I.values()) if i.state.label == "I" and i.hh_source==0])/len(inf_pop) if len(inf_pop)> 0 else np.nan
                self.stat['com_transmission_percentage']=com_transmisssion_perc
                
                #self.hh_size_comm_dists_prevalence_df=pd.DataFrame()
                self.summarystats_df2=pd.DataFrame()
                #comm transmission in age groups
                ages_com_transmitted_inf_pop=np.array([i.age for i in list(pop.I.values()) if i.state.label == "I" and i.hh_source==0])
                ages_com_transmitted_inf_pop_hist=np.histogram(ages_com_transmitted_inf_pop, bins=16, range=(0, 80))[0]
                ages_infected= np.array([i.age for i in list(pop.I.values()) if i.state.label == "I"])
                ages_infected_hist=np.histogram(ages_infected, bins=16, range=(0, 80))[0]
                #division by zero error
                #age_dist_comm_transmission=ages_com_transmitted_inf_pop_hist/ages_infected_hist
                age_dist_comm_transmission=np.array([x/y if y else np.nan for x, y in zip(ages_com_transmitted_inf_pop_hist, ages_infected_hist)])
                new_df_ages=pd.DataFrame([[t]]+[[i] for i in age_dist_comm_transmission]).transpose()
                self.age_comm_dists_prevalence_df=pd.concat([self.age_comm_dists_prevalence_df, new_df_ages], ignore_index=True)
                
                #comm transmission in hh groups
                #hh_size_comm_trans=np.array([hh_size(pop, i) for i in list(pop.I.values()) if i.state.label == "I" and i.hh_source==0])
                hh_size_comm_trans=np.array([hh_size(pop, i) for i in list(pop.I.values()) if i.state.label == "I" and i.hh_source==0])
                #display("ids of hh size comm trans inhhsize1")
                #display(np.array([[i.ID, i.groups['household']] for i in list(pop.I.values()) if i.state.label == "I" and i.hh_source==0 and hh_size(pop, i) == 1]))
                hh_size_comm_trans_hist=np.histogram(hh_size_comm_trans, bins=12, range=(1,13))[0]
                #hh_size_infected= np.array([hh_size(pop, i) for i in list(pop.I.values()) if i.state.label == "I"])
                hh_size_infected= np.array([hh_size(pop, i) for i in inf_pop])
                #display("ids of hh size infected inhhsize1")
                #display(np.array([[i.ID, i.groups['household']] for i in list(pop.I.values()) if i.state.label == "I" and hh_size(pop, i) == 1]))
                hh_size_infected_hist=np.histogram(hh_size_infected, bins=12, range=(1,13))[0]
                #division by zero error
                #hh_size_dist_comm_transmission=hh_size_comm_trans_hist/hh_size_infected_hist
                hh_size_dist_comm_transmission=np.array([x/y if y else np.nan for x, y in zip(hh_size_comm_trans_hist, hh_size_infected_hist)])
                new_df_hhs=pd.DataFrame([[t]]+[[i] for i in hh_size_dist_comm_transmission]).transpose()
               
                self.hh_size_comm_dists_prevalence_df=pd.concat([self.hh_size_comm_dists_prevalence_df, new_df_hhs],ignore_index=True)
                #display("self.hh_size_comm_dists_prevalence_df")
                #display(self.hh_size_comm_dists_prevalence_df)
                #display("new_df_hhs")
                #display(new_df_hhs)
                #display("hh_size_dist_comm_transmission")
                #display(hh_size_dist_comm_transmission)
                #display("hh_size_infected_hist")
                #display(hh_size_infected_hist)
                #display("hh_size_comm_trans_hist")
                #display(hh_size_comm_trans_hist)
                
                
                #hh with 0,1,2 cases
                number_hhs=len(pop.groups['household'])#len(pop.households) #number of total hhs in pop
                #display(list(pop.groups['household'].keys())) #list of hh ids
                #display(len(list(pop.groups['household'].keys()))) #len of hhs
                hhs_with_zerocases=0
                hhs_with_onecase=0
                hhs_with_twocases=0
                hhs_with_three_plus_cases=0
                
                for x in list(pop.groups['household'].keys()):
                    #display(x)
                    cur_hh_members= pop.groups['household'][x]
                    #display(cur_hh_members)
                    num_cases=[1 if x.state.label=="I" else 0 for x in cur_hh_members]
                    #display(num_cases)
                    if sum(num_cases)==0:
                        hhs_with_zerocases+=1
                    elif sum(num_cases)==1:
                        hhs_with_onecase+=1
                    elif sum(num_cases)==2:
                        hhs_with_twocases+=1
                    elif sum(num_cases)>=3:
                        hhs_with_three_plus_cases+=1
                          
                self.stat['hh_with_zerocase']=hhs_with_zerocases/number_hhs
                self.stat['hh_with_onecase']=hhs_with_onecase/number_hhs
                self.stat['hh_with_twocase']=hhs_with_twocases/number_hhs
                self.stat.append()
                new_stats_df=pd.DataFrame({'t':[t], 
                       'overallprevalence':[len([i for i in list(pop.I.values()) if i.state.label == "I"])/len(list(pop.I.values()))],
                       'mean_age_infected':mean_age_infected,
                        'variance_age_infected':var_age_infected,
                           'mean_hhsize_infected':mean_hhsize_infected,
                             'hh_with_zerocase':[hhs_with_zerocases/number_hhs],
                                           'hh_with_onecase':[hhs_with_onecase/number_hhs],
                                           'hh_with_twocase':[hhs_with_twocases/number_hhs],
                                           'com_transmission_percentage':com_transmisssion_perc
                    
                })
                #display(new_stats_df)
                
                self.stats_df=pd.concat([self.stats_df, new_stats_df], ignore_index=True)
                if t/self.t_per_year == round(self.et/self.t_per_year): #final saving moment #save documents into csv files
                    #display(self.fname[:-4])
                    self.stats_df.to_csv('%s_stats.csv'%(self.fname[:-4]), index=False)
                    self.age_dists_prevalence_df.to_csv('%s_agedists.csv'%(self.fname[:-4]), index=False)
                    self.hh_size_dists_prevalence_df.to_csv('%s_hhdists.csv'%(self.fname[:-4]), index=False)
                    
                    #take the summaries
                    stats_mean=self.stats_df.mean(axis=0)[1:]
                    #display(stats_mean)
                    #display(stats_mean[len(stats_mean)-1:])
                    agedist_mean=self.age_dists_prevalence_df.mean(axis=0)[1:]
                    hhdist_mean=self.hh_size_dists_prevalence_df.mean(axis=0)[1:]
                    age_comm_dist_mean=self.age_comm_dists_prevalence_df.mean(axis=0)[1:]
                    #display("self.hh_size_comm_dists_prevalence_df.mean(axis=0)")
                    #display(self.hh_size_comm_dists_prevalence_df.mean(axis=0))
                    #display("self.hh_size_comm_dists_prevalence_df.mean(axis=0)[1:]")
                    #display(self.hh_size_comm_dists_prevalence_df.mean(axis=0)[1:])
                    #display("self.hh_size_comm_dists_prevalence_df")
                    #display(self.hh_size_comm_dists_prevalence_df)
                    #display("self.hh_size_comm_dists_prevalence_df")
                    #display(self.hh_size_comm_dists_prevalence_df)
                    hh_comm_dist_mean=self.hh_size_comm_dists_prevalence_df.mean(axis=0)[1:]
                    #hh_comm_dist_mean=self.hh_size_comm_dists_prevalence_df.replace(0,np.nan).mean(axis=0)[1:]
                    #display(stats_mean)
                    #display(agedist_mean)
                    agedist_mean=agedist_mean.rename(index=lambda s: 'Age%s'%s)
                    hhdist_mean=hhdist_mean.rename(index=lambda s: 'HH%s'%s)
                    age_comm_dist_mean=age_comm_dist_mean.rename(index=lambda s: 'Age%s'%s)
                    hh_comm_dist_mean=hh_comm_dist_mean.rename(index=lambda s: 'HH%s'%s)
                    #display(hhdist_mean)
                    
                    #display(hhdist_mean.rename(index=lambda s: 'HH%s'%s))
                    self.summarystats_df=pd.concat([agedist_mean,hhdist_mean,stats_mean])
                    #display(self.summarystats_df)
                    self.summarystats_df.columns = ['name', 'value']
                    self.summarystats_df.to_csv('%s_summarystats.csv'%(self.fname[:-4]))
                    
                    #new summ stats2 csv file
                    self.summarystats_df2=pd.concat([age_comm_dist_mean,hh_comm_dist_mean,stats_mean[len(stats_mean)-1:]])
                    #display(self.summarystats_df2)
                    self.summarystats_df2.columns = ['name', 'value']
                    self.summarystats_df2.to_csv('%s_summarystats2.csv'%(self.fname[:-4]))
                    #display(self.summarystats_df)
                    
            
        # Need to flush occasionally to avoid PerformanceWarnings
        self.h5file.flush()