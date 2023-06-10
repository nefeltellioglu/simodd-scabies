
import os
from random import Random
from collections import defaultdict

import numpy as np
import pickle as pickle

from population.utils import create_path
from disease.general.sim_epi import SimEpi


def vaccinate_population_individual(sim, prop, rng):
    """
    Vaccinate (move to R) a specfied proportion of individuals.
    """
    for ind in list(sim.P.I.values()):
        if rng.random() < prop and ind.state.at_risk:
            ind.next_state = sim.disease.states['R']
            sim.disease.tick(0, ind)
            
    
def vaccinate_population_household(sim, prop, rng):
    """
    Vaccinate (move to R) individuals belonging to a specified proportion of households.
    """
    for hh in list(sim.P.groups['household'].values()):
        if rng.random() < prop and all(x.state.at_risk for x in hh):
            for ind in hh:
                ind.next_state = sim.disease.states['R']
                sim.disease.tick(0, ind)


def seed_population(sim, seed_cases):
    """
    Set the specified individual(s) to infected.
    """
    
    if not isinstance(seed_cases, list):
        seed_cases = [seed_cases]
    #sim.disease.seed_infection(0, sim.P, len(seed_cases), None, seed_cases)
    for cur_case in seed_cases:
        cur_case.next_state = sim.disease.states['I']
        sim.disease.tick(0, cur_case)
        sim.disease.I_by_age[sim.disease.cmatrix.age_map[cur_case.age]] += 1


def reset_population(sim, cases):
    """
    Reset population such that all individuals are susceptible.
    """
    
    for cur_case in cases:
        cur_case.next_state = sim.disease.states['S']
        sim.disease.tick(0, cur_case)
    sim.disease.set_counts(sim.P) 
    sim.disease.observers['cases']._reset()

    
def get_secondary_cases(sim, cur_case):
    """
    Get the secondary cases arising from a specified case.
    Leaves population unmodified (ie, reverted to susceptible)
    """
    
    seed_population(sim, cur_case)
#    sim.disease.update(1, sim.P, sim.rng)
    cases = sim.disease.check_exposure(1, sim.P, sim.rng)       
    secondary_cases = sim.disease.update_ind_states(1, sim.P)
    reset_population(sim, [cur_case] + secondary_cases)
    return secondary_cases
    

def get_primary_case(sim, cur_case, steps, rng):
    """
    A recursive function to sample a "typical infected individual" from the population.
    The returned individual is at the end of a chain of infection of length "steps".
    """
    
    if steps > 0 and cur_case:
        secondary_cases = get_secondary_cases(sim, cur_case)
        next_case = rng.sample(secondary_cases, 1)[0] if secondary_cases else None
        return get_primary_case(sim, next_case, steps-1, rng)
    else:
        return cur_case
    

def estimate_R0(p, disease_type, cmatrix, cur_seed, sim_type=SimEpi):
    """
    Estimate R0 by introducing infectious individuals randomly into 
    susceptible populations and measuring number of secondary cases.

    Also calculates proportion of household members infected by randomly
    introduced infectious individual.
    """
    
    ofile = os.path.join(p['prefix'], 'output.pickle')
    if os.path.isfile(ofile):
        f = file(ofile, 'rt')
        output = pickle.load(f)
        f.close()
        return output
    
    create_path(p['prefix'])
    p = p.copy() # create local copy of param dictionary
    disease_fname = os.path.join(p['prefix'], 'disease_%d.hd5' % cur_seed)

    rng = Random(cur_seed)
    trial_rng = np.random.RandomState(cur_seed)
    
    # setup population to have demographic structure but no disease
    p['initial_cases'] = 0
    p['external_exposure_rate'] = 0.0
    p['epi_burn_in'] = 0
    p['t_per_year'] = 364 / p['infective_mean']
    disease = disease_type(p, cmatrix, rng, disease_fname)
    sim = sim_type(p, disease, rng)
    sim.setup()
    
    secondary_case_counts = []
    hh_proportions = []
    hh_sizes = []
    #secondary_cases_by_hh_size = defaultdict(list)
    
    for i in range(p['num_trials']):
        
        # identify a suitable primary case
        primary_case = None
        while primary_case == None:
            seed_case = rng.sample(list(sim.P.I.values()), 1)[0]
            primary_case = get_primary_case(sim, seed_case, p['steps'], rng)
        
        # store household size
        hh_sizes.append(sim.P.hh_size(primary_case))
        
        # get secondary cases
        secondary_cases = get_secondary_cases(sim, primary_case)

        secondary_case_counts.append(len(secondary_cases))
        #secondary_cases_by_hh_size[sim.P.hh_size(primary_case)].append(len(secondary_cases))
        
        # record secondary infections in household
        hmates = sim.P.housemates(primary_case)
        if hmates:
            hh_cases = [x for x in hmates if x.state.label == 'I']
            hh_proportions.append(float(len(hh_cases)) / len(hmates))

    # calculate summaries
    output = {}
    output['R0_estimate'] = np.mean(secondary_case_counts)
    output['R0_variance'] = np.std(secondary_case_counts)
    output['mean_hh_proportion'] = np.mean(hh_proportions)
    output['mean_hh_size'] = np.mean(hh_sizes)
    
    f = file(ofile, 'w')
    pickle.dump(output, f)
    f.close()
    
    #for cur_size in sorted(secondary_cases_by_hh_size.keys()):
    #    cur_number = len(secondary_cases_by_hh_size[cur_size])
    #    cur_avg = np.mean(secondary_cases_by_hh_size[cur_size])
    #    print cur_size, cur_number, cur_avg    
    
    return output
