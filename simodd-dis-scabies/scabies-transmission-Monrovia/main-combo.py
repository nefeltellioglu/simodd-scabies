import sys,os
repo_path = '/Users/ntellioglu/repos/'
sys.path.append(repo_path)
sys.path.append(os.path.join(repo_path, 'simodd-pop-development'))
sys.path.append(os.path.join(repo_path, 'simodd-dis-development'))
sys.path.append(os.path.join(repo_path, 'simodd-dis-development/example'))

from population.utils import create_path

from disease.general.contact_matrix import KnownContactMatrix
from disease.general.contact_matrix import ContactMatrix
from disease.general.run import go_single
from disease.SEIR.disease_SEIR import DiseaseSIS
from disease.observers.obs_disease import DiseaseObserver
#from disease.observers.obs_age_disease import AgeDiseaseObserver #doesnt work
from disease.observers.obs_age_status import AgeStatusObserver #doesnt work
from disease.observers.obs_hh_susceptible import HHSusceptibleObserver
#from disease.observers.obs_hh_status import HHStatusObserver #doesnt work
#from disease.observers.obs_hh_disease import HHDiseaseObserver #returns error about the input formats

from disease.observers.obs_all_table import AllObserver
from disease.observers.obs_cases_table import CaseObserver
from disease.experiments.param_combo import ParamComboIt
from disease.observers.obs_summary_stats import SummaryStatsObserver
from params import p

import pandas as pd
import numpy
numpy.set_printoptions(threshold=sys.maxsize)



class DiseaseModel(DiseaseSIS):
    """
    Local version of SIR disease, adding observers and vaccines specific 
    to this set of experiments.
    """

    def __init__(self, p, cmatrix, rng, fname, mode='w'):
        super(DiseaseModel, self).__init__(p, cmatrix, rng, fname, mode)
        print(fname)
        self.add_observers(
            SummaryStatsObserver(h5file=self.h5file,fname=fname, #file name for csv documents
                     observefrom=(5*p['t_per_year'] + int(p['burn_in'] * p['t_per_year'])),#observe results from that time point
                     et=int((p['burn_in'] + p['years'][1]) * p['t_per_year']), #end time of simulation to save csv files
                     t_per_year=p['t_per_year']),)    


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# - # - MAIN  - # - #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

if __name__ == '__main__':

    # setup base parameters
    p['pop_prefix'] = p['prefix']
    p['epi_prefix'] = p['prefix']
    p['epi_burn_in'] = 0
    p['halt'] = True
    p['overwrite'] = True
    p['overwrite_cp'] = True
    p['save_cp'] = False
    p['years']= [0, 50] #simulation time
    # create_path(p['prefix'])
    # print(p)
    
    #knowncmatrix=pd.read_csv(os.path.join(repo_path, 'simodd-dis-development/example/data/Liberia_contact_matrix.csv'),header=None).to_numpy()
    #print(knowncmatrix)
    
    # (basic usage) run simulation
    # go_single(p, DiseaseModel, ContactMatrix(), p['seed'], verbose=True)

    # sweep parameters
    sweep_params = [
        {'name': 'q', 'values': [0.0001]},
        {'name': 'q_h', 'values': [0.016]},
        {'name': 'infective_mean', 'values': [90]},
    ]
    p['external_exposure_rate']= 0#5e-4
    # number of different seeds to use for each parameter combination
    p['num_runs'] = 1

    # generate parameter combinations (converting iterator to list)
    param_combos = list(ParamComboIt(p, sweep_params))
    
    # just for info, print out all prefixes (these will be used as output directories)
    for x in param_combos:
        print(x['prefix'], x['seed'])
    for i in range(len(param_combos)):
        # combo index is passed in as (only) argument
        combo_num = i#int(sys.argv[1])
        cur_params = param_combos[combo_num]
        #display(cur_params['prefix'])
        # run simulation
        #go_single(cur_params, DiseaseModel, KnownContactMatrix(givenmatrix=knowncmatrix), cur_params['seed'], verbose=True)
        go_single(cur_params, DiseaseModel, ContactMatrix(smooth=False), cur_params['seed'], verbose=False)
        