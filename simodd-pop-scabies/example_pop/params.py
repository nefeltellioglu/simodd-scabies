__author__ = 'ngeard'

p = {
    # directories
    'resource_prefix': 'data/',
    'prefix': 'output/',

    # demographic model parameters
    'fertility_parity_probs': 'fertility_age_parity_probs.dat',
    'hh_composition': 'hh_comp_lib_imputed.dat',
    'age_distribution': 'age_dist_lib.dat',
    'fertility_age_probs': 'fertility_age_probs_lib.dat',
    'death_rates_m': 'death_probs_zambia_1996_m.dat',
    'death_rates_f': 'death_probs_zambia_1996_m.dat',
    'fertility_age_rates': 'fertility_rates_lib.dat',
    
    'couple_prob': 0.08,
    'leaving_prob': 0.02,
    'divorce_prob': 0.01,
    'couple_age': 21,
    'couple_age_max': 60,
    'leaving_age': 18,
    'divorce_age': 21,
    'divorce_age_max': 60,
    'partner_age_diff': -2,
    'partner_age_sd': 2,
    'min_partner_age': 16,
    'birth_gap_mean': 270,
    'birth_gap_sd': 1,

    'pop_size': 5000,
    'growth_rate': 0.0,
    'imm_rate': 0.0,
    'immigration_prob': 0.0,

    'use_parity': False,
    'dyn_rates': False,
    'update_demog': True,
    'demo_burn': 1,
    'random_seed': False,
    'seed': 1234,
    't_per_year': 4,
    'logging': False,
    'years': 10,
    'record_interval': 1
}




