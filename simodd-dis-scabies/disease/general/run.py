import os
import inspect
import tables as tb
from random import Random
from disease.general.sim_epi import SimEpi
from population.utils import create_path


def go_single(p, disease_type, cmatrix, cur_seed, sim_type=SimEpi, verbose=False, p_string=None):
    """
    Run a single simulation (or load if previously run).

    :param p: The simulation parameters
    :type p: dict
    :param disease_type: The disease model to simulate.
    :type disease_type: :class:`DiseaseBase`
    :param cmatrix: The population contact matrix.
    :type cmatrix: :class:`ContactMatrix`
    :param cur_seed: Random seed for current experiment.
    :type cur_seed: int
    :param sim_type: Simulation class to use (current options are SimEpi or SimBDI).
    :type sim_type: :class:`SimEpi`
    :param verbose: Flag to indicate whether to write output to terminal.
    :type verbose: bool
    :param p_string: an optional string containing (e.g.) parameter values to be written to terminal (for information purposes)
    :type p_string: string

    """

    create_path(p['prefix'])
    if 'shared_start' in p:
        ss_dir = os.path.join(os.getcwd(), p.get('shared_start_dir', ''))
        #print(ss_dir)
        #print(p['shared_start'])
        #print(p['prefix'])
        #print(os.path.join(ss_dir, p['shared_start']), os.path.join(p['prefix'], os.path.basename(p['shared_start'])))
        if not os.path.exists(os.path.join(p['prefix'], os.path.basename(p['shared_start']))):
            os.symlink(os.path.join(ss_dir, p['shared_start']), os.path.join(p['prefix'], os.path.basename(p['shared_start'])))

    rng = Random(cur_seed)
    
    start_year = p['burn_in'] + p['epi_burn_in']
    year_list = [(start_year + x, start_year + y) for x, y in zip(p['years'][:-1], p['years'][1:])]

    disease_fname = os.path.join(p['prefix'], 'disease_%s_%s.hd5' % (year_list[-1]))

    # a check to remove invalid files; NB: will not remove files from partial/incomplete runs (use x for that)
    if os.path.isfile(disease_fname) and not tb.is_pytables_file(disease_fname):
        os.remove(disease_fname)

    # load disease if output file already exists, otherwise run simulation
    if os.path.isfile(disease_fname) and not p['overwrite']:
        print("@@_go_single: loading existing disease (%s%sseed=%d)..." % (
            disease_type, p_string + ";" if p_string else None, cur_seed))
        print("NB: to overwrite existing disease output, rerun with 'x' switch (eg, 'python main.py s x')")
        disease = None
        try:
            disease = disease_type(p, cmatrix, rng, disease_fname, mode='r')
            
        except tb.exceptions.NoSuchNodeError:
            # this is thrown when requested observers don't exist in disease file...
            if disease:
                disease.done(False)
            print("existing disease not complete... rerunning")
        else:
            if disease.is_complete():
                disease.done()
                return
            else:
                disease.done(False)
                print("existing disease not complete... rerunning")
    if (verbose):
        print("@@_go_single: running simulation (%sseed=%d)..." % (p_string + ";" if p_string else "", cur_seed))
    # create and set up simulation object
    disease = disease_type(p, cmatrix, rng, disease_fname, mode='w')
    sim = sim_type(p, disease, rng)
    if 'fake' not in p:
        print("here, run.py")
        sim.run(verbose)
        if (verbose):
            print("\t... simulation DONE! (%sseed=%d)..." % (p_string + ";" if p_string else "", cur_seed))
        disease.done(True)

def get_times(p, start_time=0, end_time=0, days_per_tick=50):
    """
    si and ei are start_index and end_index, i.e., start and end times in units of *timesteps*
    st and et are start_time and end_time, i.e., start and end times in units of *days*
    """
    start_year = 0
    final_year = p['years'][-1] if end_time is 0 else p['years'][-2] + end_time / 364.0
    #    print "final_year", final_year
    end_year = final_year - p['years'][-2]
    times = {
        'si': int(p['t_per_year'] * start_year),
        'ei': int(p['t_per_year'] * end_year),
        'st': int((p['burn_in'] + p['epi_burn_in'] + p['years'][-2]) * p['t_per_year']) - start_time,
        'et': int((p['burn_in'] + p['epi_burn_in'] + final_year) * p['t_per_year']) - start_time,
        't_per_year': p['t_per_year'],
        'years_per_tick': 1,
        'days_per_tick': days_per_tick
    }
    #    print p['burn_in'], p['epi_burn_in'], p['years']
    #    print times
    return times


def output_single(p, disease):
    """
    Write all output for a single simulation.

    :param p: The simulation parameters.
    :type p: dict
    :param disease: The completed disease simulation object.
    :type disease: :class:`DiseaseBase`

    """
    print("Processing data and generating output plots...")
    times = get_times(p)  # , end_time=75, days_per_tick=15)
    create_path(p['prefix'])

    # Try calling output all for each observer
    for cur_label, cur_observer in disease.observers.items():
        if inspect.ismethod(cur_observer.output_all):
            print("Processing %s output" % cur_label)
            cur_observer.output_all(p, times)
        else:
            print("Observer %s doesn't implement output_all!" % cur_label)


