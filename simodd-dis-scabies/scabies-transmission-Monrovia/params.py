__author__ = 'ngeard'

p = {
    # directories
    'resource_prefix': 'data/',
    'prefix': 'output/',

    # demographic model parameters
    'fertility_age_rates': 'fertility_rates.dat',
    'hh_composition': 'hh_comp_zambia_2000.dat',
    'age_distribution': 'age_dist_lib.dat',
    'fertility_parity_probs': 'fertility_age_parity_probs.dat',
    'fertility_age_probs': 'fertility_age_probs.dat',
    'death_rates_m': 'death_probs_zambia_1996_m.dat',
    'death_rates_f': 'death_probs_zambia_1996_m.dat',
    'birth_num_file': 'birth_num.dat',
    
    'couple_prob': 0.039, #based on zambia sample population:0.08
    'leaving_prob':0.008,#0.00005, #zambia sample population:0.005
    'divorce_prob':0.001, #based on zambia sample population: 0.001
    'couple_age':15, #is this min couple age? #based on zambia sample population
    'couple_age_max': 60,
    'leaving_age': 18,
    'divorce_age': 18, #is this min divorcing age? #based on zambia sample population
    'divorce_age_max': 60, #based on zambia sample population
    'partner_age_diff':2, #should this be -2 or 2? #based on zambia sample population
    'partner_age_sd': 2, #based on zambia sample population
    'min_partner_age': 15, #based on zambia sample population
    'birth_gap_sd': 0, #based on zambia sample population
    'birth_gap_mean':270, #based on zambia sample population
    

    'pop_size': 5000,
    'growth_rate': 0.01,
    'imm_rate': 0.0,
    'immigration_prob': 0.0,

    'preg': False,
    'use_parity': False,
    'dyn_rates': False,
    'update_demog': True,

    # contact model parameters
    'cm_gauss': False,
    'cm_smooth': False,
    'sigma_2': 10.0,
    'epsilon': 0.8,
    'phi': 0.7,
    'cm_update_years': 5,

    # disease model parameters (durations in days)
    'infective_mean': 7 * 1,
    'norecovery_cases_frac': 0.1,
    'alpha': 1.0,
    'sf_amp': 0.0, #seasonal 
    'q': 0.00,
    'q_h': 0.3,
    'k': 3,
    'external_exposure_rate': 5e-6,

    'random_seed': False,
    'seed': 1234,
    't_per_year': 52,
    'years': [0, 10],
    'burn_in': 100,
    'burn_in_t_per_year': 4,
    'epi_burn_in': 100,

    'halt': False,

    # run parameters
    'num_runs': 1,
    'initial_cases': round(0.093 * 13635), #lib data prevalence is 9.3%
    'output_list': ['all'],
    'save_cp': True,
    'logging': False,
}

p['demo_burn'] = p['burn_in'] + p['epi_burn_in']



