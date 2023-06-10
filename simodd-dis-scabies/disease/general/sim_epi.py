"""
sim.py

Simulation loop for epidemic simulations.

Eventually this should handle population updating via a "pop_updater"
class analogous to the "disease_updater" class, and stop inheriting pop.simulation.

TODO: the firstborn/subsequent counts are experiment-specific and should be
handled via an observer, rather than as currently.  Working out how best to do this
may be tricky as it requires extending observer capability to demographic updating.
"""
import time
import os
import pickle as pickle
import numpy as np

from population.pop_io import load_pop, save_pop
from population.simulation import Simulation
from disease.general.ind_epi import IndEpi


class SimEpi(Simulation):
    def __init__(self, p, disease, rng, ind_type=IndEpi):
        super(SimEpi, self).__init__(p, ind_type, h5file=disease.h5file, create_pop=False)
        self.disease = disease
        self.rng = rng
        self.disease.firstborn_count = 0
        self.disease.subsequent_count = 0

    def _init_contact_matrix(self, t=0):
        # initialise contact matrix
        # print "Initialising contact matrix; t=%d" % t
        # print self.P.group_size_dist('household', max_size=20, norm=False)
        # print self.P.get_mean_housemates_by_age(cutoffs=[0, 1, 5, 15, 19, 59])
        self.disease.cmatrix.init_age_classes(self.P)
        # self.disease.cmatrix.init_a_levels(self.P)
        if self.p['cm_gauss']:
            self.disease.cmatrix.init_contact_matrix_gaussian(
                self.p['epsilon'], self.p['sigma_2'])
        else:
            self.disease.cmatrix.init_contact_matrix(self.p['epsilon'])
        #print(repr(self.disease.cmatrix.C))
        #print(self.disease.cmatrix.EC)

    def _init_contact_network(self):
        net_filename = os.path.join(self.p['pop_prefix'],
                                    'cnetwork_%d.pop' % self.p['pop_size'])

        if os.path.isfile(net_filename):
            print("Loading previously created contact network")
            self.disease.cnetwork = pickle.load(open(net_filename, 'rb'))
            print("Loaded contact network")
        else:
            print("No contact network found in %s; creating new contact network..."
                  % self.p['pop_prefix'])

            self.disease.cnetwork.create_network(self.P, self.disease.cmatrix, self.rng)
            pickle.dump(self.disease.cnetwork, open(net_filename, 'wb'))

    def _demo_burn_in(self):
        """
        Load burnt in population if exists.

        Run the demographic component of the simulation for specified number
        of years (typically 100) using (typically!) annual updating steps 
        in order to generate a starting population that avoids the artifactual
        traits of the bootstrap population
        """

        # if a specific location for the burnt in population isn't 
        # specified, look for (and save) it in the output directory
        if 'pop_prefix' not in self.p:
            print("WARNING!: pop_prefix not specified!")
            self.p['pop_prefix'] = self.p['prefix']

        pop_filename = os.path.join(self.p['pop_prefix'],
                                    'burn_in_%d.pop' % self.p['pop_size'])

        if os.path.isfile(pop_filename):
            #print("Loading previously burnt in population")
            self.P = load_pop(pop_filename)
            #print("Loaded (%d individuals in %d households)" % (len(self.P.I), len(self.P.groups['household'])))
        else:
            print("No burn in population found in %s; creating new population..."
                  % self.p['pop_prefix'])

            timesteps = self.p['burn_in'] * (self.p['burn_in_t_per_year'])
            real_t_per_year = self.p['t_per_year']
            self.p['t_per_year'] = self.p['burn_in_t_per_year']
#            self.load_demographic_data()
            success = False
            while not success:
                self.create_population()
                for i in range(timesteps):
                    # print "burn in; year", i
                    b, d, i, b2 = self.update_all_demo(i)
                    if b == "error":
                        success = False
                        break
                else:
                    success = True
            self.p['t_per_year'] = real_t_per_year
            self._load_demographic_data()
            print("Burn in complete; saving population...")
            save_pop(self.P, pop_filename)
            print("Population saved.")

    def _epi_burn_in(self, start_year, verbose=True):
        """
        Load burnt in population if exists.

        Otherwise, run an epidemic for specified number of years (i.e., 100) to 
        produce a stable population states.  Save the resulting population.
        """

        # if a specific location for the burnt in population isn't 
        # specified, look for (and save) it in the output directory
        if 'epi_prefix' not in self.p:
            print("WARNING!: epi_prefix not specified!")
            self.p['epi_prefix'] = self.p['prefix']
        #May42021 Nefel-> changed name of epiburnin file
        #epi_filename = os.path.join(self.p['epi_prefix'],
        #                            'burn_in_epi_%d.pop' % self.p['pop_size'])
        #if not self.p['external_exposure_rate']:
        epi_filename = os.path.join(self.p['epi_prefix'],
                                    'burn_in_epi_%d_q%s_qh%s_inf%s.pop' %(self.p['pop_size'],self.p['q'],self.p['q_h'],self.p['infective_mean']) )
        #else:
         #   epi_filename = os.path.join(self.p['epi_prefix'],
         #                           'burn_in_epi_%d_q%s_qh%s_inf%s_exposure%s.pop' %(self.p['pop_size'],self.p['q'],self.p['q_h'],self.p['infective_mean'], self.p['external_exposure_rate']) )
        if os.path.isfile(epi_filename):
            print("Loading previously burnt in epidemic")
            self.P = load_pop(epi_filename)
            self.disease.set_counts(self.P)
        else:
            print("No burn in epidemic found in %s; creating new epidemic."
                  % self.p['epi_prefix'])
            self._demo_burn_in()
            print("Demographic burn in complete/loaded; running epidemic burn in...")

            self._init_contact_matrix()
            # seed infection, switch observers off, burn in and save population
            self.disease.seed_infection(start_year * self.p['t_per_year'], self.P,
                                        self.p['initial_cases'], self.rng)
            self.disease.obs_on = False
            self._main_loop(self.p['burn_in'], self.p['epi_burn_in'], verbose)
            self.disease.obs_on = True
            print("Epidemic burn in complete; saving population...")
            save_pop(self.P, epi_filename)
            print("Population saved.")

    def setup(self, start_year=0, verbose=True, seed_inds=None):
        """
        Setup population for start of simulation.
        
        If no epidemic burn in is required, run (or load) demographic burn in  
        and seed infection.  Otherwise, run (or load) epidemic burn in.
        """

        if self.p['epi_burn_in'] <= 0:
            self._demo_burn_in()
            self.disease.seed_infection(start_year * self.p['t_per_year'], self.P,
                                        self.p['initial_cases'], self.rng, seed_inds)
        else:
            self._epi_burn_in(start_year, verbose)

        # if run procedure is *always* called via run_cp, then this can be removed (?)
        if self.disease.cnetwork:
            # (TODO: this is a bit ugly)
            print("Initialising contact matrix (setup)...")
            self._init_contact_matrix(start_year * self.p['t_per_year'])
            print("Initialising contact network...")
            self._init_contact_network()

    def _main_loop(self, year_begin, years, verbose=False):
        """
        Run simulation.
        """
        t_begin = int(year_begin * self.p['t_per_year'])
        t_end = int((year_begin + years) * self.p['t_per_year'])

        if verbose:
            self.start_time = time.time()
            self.print_column_labels()
            self.print_pop_numbers(t_begin)

        self.disease.update_observers(t_begin, disease=self.disease, pop=self.P,
                                      cases=[],
                                      boosting=[],
                                      introduction=False,
                                      new_I=[], rng=self.rng)

        for t in range(t_begin + 1, t_end + 1):
            # update demography (if required)
            if self.p['update_demog']:  # and t%52==0:
                births, deaths, imms, birthdays = self.update_all_demo(t)  # *52)
                firstborns = len([x for x in births if x.birth_order == 1])
                self.disease.firstborn_count += firstborns
                self.disease.subsequent_count += (len(births) - firstborns)
                # oh wow, really don't need to be doing THIS all the time!
                # a) no point unless also reinitialising contact matrix;
                # b) no point at all unless population structure is changing over time
                # self.disease.cmatrix.update_age_classes(
                #        births, deaths, imms, birthdays)
                # if self.p['dyn_rates'] and t % (self.p['cm_update_years'] * self.p['t_per_year']) == 0:
                if self.p['cm_update_years'] > 0 and t % (self.p['cm_update_years'] * self.p['t_per_year']) == 0:
                    self._init_contact_matrix(t)
                self.disease.bd_update(t, births, deaths, imms, self.rng)
            
            
            # update disease
            if self.disease.update(t, self.P, self.rng):
                if verbose:
                    self.print_pop_numbers(t)
                break  # update returns true if halting upon fade out

            if verbose:
                self.print_pop_numbers(t)

        if verbose:
            self.print_column_labels()
            self.end_time = time.time()
            print("time:", self.end_time - self.start_time)

    def run_cp(self, start_year, years, observers, verbose, prefix,final_year):
        """
        A version of run (above) that checks for the existence of a
        checkpoint file before running, and loads that by preference.
        Otherwise, run is called normally, and a checkpoint file is 
        saved (in directory specified by prefix).
        
        'observers' is a Boolean value specifying whether observers
        should be switched on or off for this period.
        """
        #nefel
        #print("run_cp is called")
        cp_filename = os.path.join(prefix, 'cp_%d-%d.pop' % (start_year, start_year + years))
        #print("Current checkpoint: %s (observers=%s)" % (cp_filename, observers))
        if os.path.isfile(cp_filename) and not self.p['overwrite_cp']:
            #print("Loading previous checkpoint (years %d-%d)" % (start_year, start_year + years))
            self.P = load_pop(cp_filename)
            self.disease.set_counts(self.P)
            start_year = start_year + years
            years = final_year - start_year
            #print("%s \t %s \t %s"%(start_year,final_year,years ))
        else:
            print("No checkpoint found for year %d; running..." % start_year)

        self._init_contact_matrix(start_year * self.p['t_per_year'])

        # switch observers off, burn in and save population
        self.disease.obs_on = observers
        self.disease.firstborn_count = 0
        self.disease.subsequent_count = 0
        self._main_loop(start_year, years, verbose)
        if self.p['save_cp']:
            cp_filename = os.path.join(prefix, 'cp_%d-%d.pop' % (start_year, start_year + years))
            save_pop(self.P, cp_filename)

    def run(self, verbose):
        #print("run is called")
        start_year = self.p['burn_in'] + self.p['epi_burn_in']
        #year_list = [(start_year + x, start_year + y) for x, y in zip(self.p['years'][:-1], self.p['years'][1:])]
        #nefel
        year_list = [(start_year + self.p['years'][0], start_year + self.p['years'][0] + y) \
           for y in range(1,self.p['years'][1] - self.p['years'][0]+1)] #[(start_year , start_year + self.p['years'][0])] + \
        start_index = -1
        # starting at end, test for existence of each checkpoint, to identify where we need to start from
        #print("start index: %s"%start_index)
        for i, cur_years in enumerate(reversed(year_list)):
            cp_filename = os.path.join(self.p['prefix'], 'cp_%d-%d.pop' % (cur_years[0], cur_years[1]))
            #print("Testing for checkpoint: %s" % cp_filename)
            #display("Testing for checkpoint: %s" % cp_filename)
            
            # if this checkpoint exists, break and run forward from here, otherwise check for previous...
            if os.path.isfile(cp_filename) and not self.p['overwrite']:
                start_index = len(year_list) - i
                print(i)
                break

        # if final checkpoint exists, our work here is done...
        # this implies that disease should exist, so should never occur in practice
        if start_index == len(year_list):
            return

        # if none of our checkpoints existed, run setup to make sure populations exist
        if start_index == -1:
            self.setup(verbose=verbose)
        #    start_index = 1
        #print("start_index %s"%start_index)
        # load/run each checkpoint in turn (starting with latest existing one to get population loaded... need fixing!)
        if start_index == -1:
            cur_years =  cur_years = [(start_year + x, start_year + y) for x, y in zip(self.p['years'][:-1], self.p['years'][1:])][0]
        else: 
            cur_years = year_list[(start_index - 1): start_index][0]
            display(cur_years)
        #for cur_years in year_list[start_index - 1: start_index]:
            # for efficiency, only turn on observers for final period
        observers = True #(cur_years == year_list[-1])
        cp_prefix = self.p['prefix']
        #print("years %s_%s"%(cur_years[0], cur_years[1]))
        final_year = self.p['burn_in'] + self.p['epi_burn_in'] + self.p['years'][1]
        self.run_cp(cur_years[0], cur_years[1] - cur_years[0], observers, verbose, cp_prefix,final_year)

    # OUTPUT AND HELPER FUNCTIONS ###########################################

    def print_column_labels(self):
        print(self.disease.state_labels())
        print(''.join(['%7s' % x
                       for x in
                       (['t', 't(y)', 't(d)'] + self.disease.state_labels())]))

    def print_pop_numbers(self, t):
        print(''.join(['%7d' % x
                       for x in
                       ([t, t // self.p['t_per_year'],
                         t % self.p['t_per_year'] * (364.0 / self.p['t_per_year'])]
                        + self.disease.state_counts())]))
