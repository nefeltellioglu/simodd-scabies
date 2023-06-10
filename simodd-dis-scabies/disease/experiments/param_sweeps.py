import inspect
import os
from random import Random
from itertools import product
from mpi4py import MPI
import tables as tb

import numpy as np
import pickle as pickle
from disease.general.contact_matrix import ContactMatrix
from disease.general.sim_epi import SimEpi
from population.utils import create_path

# An iterator object for parallelising across multiple parameter combinations


class ParamComboIt(object):
    """
    An iterator object for parameter sweeps.  Generates dictionaries containing
    unique combinations of parameter values and random seeds.  Each dictionary
    stores sufficient information to be farmed out to independent processes
    using :func:`go_sweep_mpi`.

    :param p: The simulation parameters
    :type p: dict
    :param p_values: A list of dictionaries containing parameters and values to sweep over
    :type p_values: list
    :param disease_type: The disease model to simulate.
    :type disease_type: :class:`DiseaseBase`
    :param cmatrix: The population contact matrix
    :type cmatrix: :class:`ContactMatrix`
    :param use_seeds: Flag to specify whether iterator should include seeds (set to False if wanting to aggregate over seeds)
    :type use_seeds: bool
    :param output_list: A list of output type flags to aggregate
    :type output_list: list

    """

    def __init__(self, p, p_values, disease_type=None, cmatrix=None, use_seeds=True):
        self.p = p
        self.disease_type = disease_type
        self.cmatrix = cmatrix
        self.output_list = p['output_list'] if 'output_list' in p else ['all']
        # create a RNG for seeds (we will use same seeds for each param combo)
        # the use_seeds flag allows seeds to be ignored
        if use_seeds:
            seed_rng = Random() if p['random_seed'] else Random(p['seed'])
            self.seeds = [seed_rng.randint(0, 99999999) for _ in range(p['num_runs'])]
        else:
            self.seeds = [p['seed']]
        self.base_prefix = p['prefix']
        self.p_values = p_values
        # generate all necessary parameter combinations (indexed by sample_ID)
        self.combos = [a for a in product(*[list(range(len(x['values']))) for x in p_values])]

        self.combo_index = 0
        self.seed_index = 0

    def print_combos(self):
        """
        Print all parameter combinations.

        """

        for cur_combo in self.combos:
            print(("".join(['{0:>10}'.format(self.p_values[i]['values'][x])
                           for i, x in enumerate(cur_combo)])))

    def __iter__(self):
        return self

    def __next__(self):
        """
        Get the next parameter combination object.

        :returns:parameter combination object (dictionary).
        """

        if self.combo_index == len(self.combos):
            raise StopIteration
        else:
            #print "combo:", self.combo_index + 1, ", (", len(self.combos), ")"
            #print "seed:", self.seed_index + 1, ", (", len(self.seeds), ")"

            # make a local copy of params to modify
            p = self.p.copy()
            # set parameter values for this combination
            cur_combo = self.combos[self.combo_index]
            output_dir = "params_"
            param_strings = []
            for i, cur_p_name in enumerate([x['name'] for x in self.p_values]):
                p[cur_p_name] = self.p_values[i]['values'][cur_combo[i]]
                param_strings.append("%s=%s" % (cur_p_name,
                                                self.p_values[i]['values'][cur_combo[i]]))
            output_dir += "_".join(param_strings)
            p['prefix'] = os.path.join(self.base_prefix, output_dir)
            # only use seed directories if running more than one seed
            if len(self.seeds) > 1:
                p['prefix'] = os.path.join(p['prefix'], 'seed_%d' % self.seed_index)

            cur_seed = self.seeds[self.seed_index]
            #print "cur_seed", cur_seed

            if self.seed_index < len(self.seeds) - 1:
                #print "incrementing seed"
                self.seed_index += 1
            else:
                if self.combo_index < len(self.combos):
                    #print "incrementing combo"
                    self.combo_index += 1
                    #print "resetting seed"
                    self.seed_index = 0

            return {
                'disease_type': self.disease_type,
                'cmatrix': self.cmatrix,
                'seed_index': self.seed_index,
                'seed': cur_seed,
                'p': p,
                'p_string': ";".join(param_strings),
                'output_list': self.output_list
            }


class SweepResults(object):
    """
    Stores a multidimensional array mapping parameter combinations to result indices.
    The aim is to make it straightforward to pull out one- and two-dimensional parameter
    maps.
    """

    def __init__(self, sweep_values, combos, num_seeds):
        self.sweep_values = sweep_values
        self.num_seeds = num_seeds
        # create a mapping from a parameter hypercube to an output
        self.indices = np.zeros(tuple([len(x['values']) for x in sweep_values]), dtype=int)
        # populate array with indices
        for i, cc in enumerate(combos):
            pt = "self.indices" + "".join(["[cc[%d]]" % x for x in range(len(cc))]) + "=%d" % i
            exec (pt)
        # create a mapping from parameter name to axis index
        self.index_map = dict([(x['name'], i) for i, x in enumerate(sweep_values)])

    def retrieve_results(self, fixed=None):
        """
        retrieve results:          
        - provide a list of (param, value) pairs.
        - use index_map and sweep_values to get (param_index, value_index) pairs
        - build a string representing array access (assume ':' for any non-specified params)
        - execute and return result subset

        """

        if not fixed:
            fixed = {}
        subset = None
        pt = "subset = self.indices["
        pt += ",".join([
            ":" if x['name'] not in fixed else "%d" % x['values'].index(fixed[x['name']])
            for x in self.sweep_values])
        pt += "]"
        exec (pt)
        # if sweep is 1D, need to force return value to be a list (of length one)
        return [subset] if isinstance(subset, int) else subset.tolist()

    def retrieve_results_sweep(self, sweep):
        """
        Return a list of results obtained by iteratively holding each value of the sweep param fixed,
        and sweeping over the remaining parameters.
        NB: for a param hypercube of n dimensions, will return slices of dimension n-1
        eg, for a 2D set of param combinations, will return a list of results;
        for a 3D set of param combinations, will return a matrix of results, etc.
            
        TODO: should return values as well as results (as per below)
        """
        return [self.retrieve_results({sweep: x})
                for x in self.sweep_values[self.index_map[sweep]]['values']]

    def retrieve_results_multi_sweep(self, cur_param):
        """
        Return a list of results obtained by holding sweeping over cur_p for all combinations of other parameters.
        NB: will *always* return a *list* of results (one per parameter combination).
        """
        # generate combinations of remaining sweep parameters
        other_names = [x['name'] for x in self.sweep_values if x['name'] != cur_param]
        other_combos = [a for a in product(*[list(range(len(x['values'])))
                                             for x in self.sweep_values if x['name'] != cur_param])]

        values = []
        results = []

        # for each combination of other sweep parameters
        for cur_combo in other_combos:
            # build a dictionary of current values
            cur_values = {}
            for i, cur_other in enumerate(other_names):
                cur_values[cur_other] = self.sweep_values[self.index_map[cur_other]]['values'][cur_combo[i]]
            values.append(list(cur_values.values()))
            results.append(self.retrieve_results(cur_values))

        return values, results


def go_sweep_mpi(p, sweep_values, disease_type, cmatrix, process_fn):
    """
    Distribute a series of simulations over multiple processors (using MPI).

    :param p: The simulation parameters
    :type p: dict
    :param sweep_values: A list of dictionaries containing parameters and values to sweep over
    :type sweep_values: list
    :param disease_type: The disease model to simulate.
    :type disease_type: :class:`DiseaseBase`
    :param cmatrix: The population contact matrix
    :type cmatrix: :class:`ContactMatrix`
    :param process_fn: The function to run on each processor.
    :type process_fn: function pointer

    """

    # produce a generator object for all possible parameter/seed combinations
    param_combo_it = ParamComboIt(p, sweep_values, disease_type, cmatrix)
    sweep_results = SweepResults(sweep_values, param_combo_it.combos, p['num_runs'])

    # show parameter combinations
    print(("@@_go_sweep: About to sweep over the following parameter", \
        "combinations using function %s:" % process_fn.__name__))
    print(("".join(["%10s" % x['name'] for x in sweep_values])))
    param_combo_it.print_combos()

    comm = MPI.COMM_WORLD
    mode = MPI.MODE_WRONLY|MPI.MODE_CREATE
    fh = MPI.File.Open(comm, 'logfile.log', mode)
    fh.Set_atomicity(True)
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        params_all = [x for x in param_combo_it]
        params_slices = [[] for _ in range(size)]
        for i, x in enumerate(params_all):
            # noinspection PyTypeChecker
            params_slices[i % size].append(x)
    else:
        params_slices = None

    params_local = comm.scatter(params_slices, root=0)
    results_local = [process_fn(cur_params, fh, rank) for cur_params in params_local]
    results_combined = comm.gather(results_local, root=0)

    fh.Sync()
    fh.Close()

    if rank == 0:
        output = [x for y in results_combined for x in y]
        return sweep_results, output
    else:
        return None, None


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
        print(ss_dir)
        print((p['shared_start']))
        print((p['prefix']))
        print((os.path.join(ss_dir, p['shared_start']), os.path.join(p['prefix'], os.path.basename(p['shared_start']))))
        if not os.path.exists(os.path.join(p['prefix'], os.path.basename(p['shared_start']))):
            os.symlink(os.path.join(ss_dir, p['shared_start']), os.path.join(p['prefix'], os.path.basename(p['shared_start'])))

    rng = Random(cur_seed)

    start_year = p['burn_in'] + p['epi_burn_in']
    year_list = [(start_year + x, start_year + y) for x, y in zip(p['years'][:-1], p['years'][1:])]

    disease_fname = os.path.join(p['prefix'], 'disease_%s_%s.hd5' % (year_list[-1]))

    # a check to remove invalid files; NB: will not remove files from partial/incomplete runs (use x for that)
    if os.path.isfile(disease_fname) and not tb.isPyTablesFile(disease_fname):
        os.remove(disease_fname)

    # load disease if output file already exists, otherwise run simulation
    if os.path.isfile(disease_fname) and not p['overwrite']:
        print(("@@_go_single: loading existing disease (%s%sseed=%d)..." % (
            disease_type, p_string + ";" if p_string else None, cur_seed)))
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

    print(("@@_go_single: running simulation (%sseed=%d)..." % (p_string + ";" if p_string else "", cur_seed)))
    # create and set up simulation object
    disease = disease_type(p, cmatrix, rng, disease_fname, mode='w')
    sim = sim_type(p, disease, rng)
    if 'fake' not in p:
        sim.run(verbose)
        print(("\t... simulation DONE! (%sseed=%d)..." % (p_string + ";" if p_string else "", cur_seed)))
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
    for cur_label, cur_observer in list(disease.observers.items()):
        if inspect.ismethod(cur_observer.output_all):
            print(("Processing %s output" % cur_label))
            cur_observer.output_all(p, times)
        else:
            print(("Observer %s doesn't implement output_all!" % cur_label))


def process_diseases_parallel(obj, fh, rank):
    """
    Wrapper function for :func:`go_single` when conducting parameter sweeps.

    :param obj: The parameter combination object returned by :class:`ParamComboIt`
    :type obj: dict
    """
    print("@@@_process_diseases_parallel")
    msg = '[%d] %s -- %d (%d) :: begin\n' % (rank, obj['p_string'], obj['seed_index'], obj['seed'])
    fh.Write_shared(msg)
    go_single(obj['p'], obj['disease_type'], obj['cmatrix'], obj['seed'], p_string=obj['p_string'])
    msg = '[%d] %s -- %d (%d) :: end\n' % (rank, obj['p_string'], obj['seed_index'], obj['seed'])
    fh.Write_shared(msg)


def output_diseases_parallel(obj):
    """
    Wrapper function for :func:`output_single` when conducting parameter sweeps.

    :param obj: The parameter combination object returned by :class:`ParamComboIt`
    :type obj: dict
    """

    print("@@@_output_diseases_parallel")
    p = obj['p']
    cur_seed = obj['seed']
    rng = Random(cur_seed)
    disease_fname = os.path.join(p['prefix'], 'disease.hd5')
    print(("opening disease file:", disease_fname))
    disease = obj['disease_type'](p, obj['cmatrix'], rng, disease_fname, mode='r')
    output_single(p, disease)
    disease.done()


def go_output_aggregate(obj):
    """
    Function for extracting numerical outputs from observers, suitable to be aggregated
    across multiple seeds.

    The output_list specifies which outputs to generate (may wish to switch some off for
    efficiency); possible values are:
    - inc: general incidence
    - age: age-specific incidence
    - hh: household size-specific incidence
    - first: age of first infection
    - hh_frac: fraction of infections attributable to hh transmission
    - susc: presence/size of susceptibility clusters
    - hh_imm: proportion of hh members immune

    :param obj: The parameter combination object returned by :class:`ParamComboIt`
    :type obj: dict
    :return: Dictionary of extracted values.

    """
    print("@@@_go_output_aggregate")
    disease = obj['disease_type'](obj['p'], obj['cmatrix'], Random(obj['seed']),
                                  os.path.join(obj['p']['prefix'], 'disease.hd5'), mode='r')
    if not disease.is_complete():
        print("--- disease file not complete... skipping!")
        return

    p = obj['p']
    times = get_times(p)
    output_list = obj['output_list']
    delta_t = 5

    output = {}
    for cur_output_type in output_list:
        if cur_output_type in ['inc', 'all']:
            # store annual incidence
            output['annual_inc'] = disease.observers['cases'].get_annual_incidence(
                times, times['t_per_year'] * delta_t)

            # store annual incidence -- community transmission only
            output['annual_inc_comm'] = disease.observers['cases'].get_annual_incidence(
                times, times['t_per_year'] * delta_t, 'comm')

        if cur_output_type in ['age', 'all']:
            # store incidence by age snapshots (binned)
            #    age_bins = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 50]
            age_bins = [0, 1, 5, 10, 20, 40, 60]
            incs_snapshots = disease.observers['cases'].get_incidence_by_age_bin(
                age_bins + [100], times['st'], times['et'], times['t_per_year'] * 20)

            for i, cur_snapshot in enumerate(incs_snapshots['I']):
                #print type(cur_snapshot)
                output['inc_by_age_bin_ss_%d' % i] = list(cur_snapshot)

            # store incidence by age snapshots (non-binned)
            incs_snapshots = disease.observers['cases'].get_incidence_by_age_snapshots(
                times['st'], times['et'], 5)

            for i, cur_snapshot in enumerate(incs_snapshots['I']):
                output['inc_by_age_ss_%d' % i] = list(cur_snapshot)

        if cur_output_type in ['hh', 'all']:
            incs_snapshots = disease.observers['cases'].get_incidence_by_hh_size_snapshots(
                times['st'], times['et'], 5, max_hh_size=10)

            for i, cur_snapshot in enumerate(incs_snapshots):
                output['inc_by_hh_size_ss_%d' % i] = list(cur_snapshot)

        if cur_output_type in ['first', 'all']:
            # store average age at first infection
            avg_age, avg_age_err = list(zip(*disease.observers['cases'].get_series(
                disease.observers['cases'].get_avg_age_first_infection,
                times['st'], times['et'], times['t_per_year'] * delta_t)))
            output['avg_age_first_infection'] = list(avg_age)

        if cur_output_type in ['hh_frac', 'all']:
            # store household fraction snapshots
            hh_fracs, hh_fracs_err = list(zip(*disease.observers['cases'].get_hh_frac_series(
                times['st'], times['et'], times['t_per_year'] * delta_t)))
            output['hh_frac'] = list(hh_fracs)

            # store household fraction by household size
            hh_hh_frac_values = disease.observers['cases'].get_hh_hh_frac_series(
                times['st'], times['et'], times['t_per_year'] * delta_t)
            for x in hh_hh_frac_values:
                output['hh_hh_frac_%d' % x] = hh_hh_frac_values[x]

        if cur_output_type in ['susc', 'all']:
            # get household susceptible cluster sizes
            hh_at_risk = disease.observers['cases'].get_hh_at_risk_clustering(
                times['st'], times['et'], times['t_per_year'] * delta_t)
            series_by_hh_size = list(zip(*hh_at_risk))[2:]  # transform to be time series, skip hh of size 0 and 1

            for i, cur_size in enumerate(series_by_hh_size):
                output['hh_susc_cluster_%d' % i] = cur_size

        if cur_output_type in ['hh_imm', 'all']:
            # store proportion immune by household size
            hh_susc_props = disease.observers['hh_immunity'].get_pos_props(window=5 * 4)

            for i, cur_size in enumerate(hh_susc_props):
                output['hh_immunity_%d' % i] = cur_size

    disease.done()
    return output


def aggregate_over_seeds(output, num_seeds, lims=None):
    """
    Aggregate output values from individual simulations (using identical parameter values)
    across multiple random seeds.

    For each output that refers to an integer or float output value,
    the mean and standard deviation over runs are calculated.
    
    Should now handle lists (of ints and floats at least) as well, calculating mean/SD over the appropriate axis.

    :param output: A list of dictionaries of output values, keyed by output name
    :type output: list
    :param num_seeds: The number of random seeds used for each parameter combination.
    :type num_seeds: int
    :param lims: can be used to only include values within the specified range (ie to exclude error values)
    :type lims: tuple
    :returns: A list of dictionaries of aggregated values (one per parameter combination)
    
    """

    agg_output = []
    for i in range(0, len(output), num_seeds):
        # remove output entries for this param combo where there is no output
        # (ie due to a simulation that didn't complete)
        cur_output_src = [x for x in output[i:i+num_seeds] if x is not None]
        cur_agg_output = {}
        for k in list(cur_output_src[0].keys()):
            if isinstance(cur_output_src[0][k], (int, float, list, tuple, np.ndarray)):
                if lims:
                    cur_values = [x[k] for x in cur_output_src if lims[0] <= x[k] <= lims[1]]
                else:
                    cur_values = [x[k] for x in cur_output_src]

                print(cur_values)
                cur_agg_output['mean_%s' % k] = np.mean(cur_values, axis=0)
                cur_agg_output['sd_%s' % k] = np.std(cur_values, axis=0)
        agg_output.append(cur_agg_output)
    return agg_output


def concatenate_over_seeds(output, num_seeds):
    """
    Similar to :func:`aggregate_over_seeds`, except that output values are concatenated,
    rather than having mean/SD calculated

    NB: as no numeric processing is done, the types of outputs are irrelevant.

    :param output: A list of dictionaries of output values, keyed by output name
    :type output: list
    :param num_seeds: The number of random seeds used for each parameter combination.
    :type num_seeds: int
    :returns: A list of dictionaries of concatenated values (one per parameter combination)

    """
    concat_output = []
    for i in range(0, len(output), num_seeds):
        cur_output = {}
        for k in list(output[i].keys()):
            cur_output[k] = [x[k] for x in output[i:i + num_seeds]]
        concat_output.append(cur_output)
    return concat_output


def agg_output_pickle(p, sweep_params, dis_class, file_prefix='agg_', path_prefix='.'):
    """
    Calculate aggregated output values, and then cache (using pickle) the result.  This can speed things up
    when developing plots (as potentially expensive observer output analyses don't need to be
    continually rerun).

    :param p: The simulation parameters
    :type p: dict
    :param sweep_params: A list of dictionaries containing parameters and values to sweep over
    :type sweep_params: list
    :param dis_class: The disease model to simulate.
    :type dis_class: :class:`DiseaseBase`
    :param file_prefix: File prefix to use for pickled files.
    :type file_prefix: string
    :param path_prefix: Path prefix to use for pickled files.
    :type path_prefix: string
    :returns: A list of dictionaries of aggregated values (one per parameter combination)

    """

    fname = file_prefix
    fname += '_'.join(['%s=%g' % (x['name'], x['values'][0]) for x in sweep_params])
    fname += '.dat'

    aggfile = os.path.join(path_prefix, fname)

    print(("checking for existence of %s" % aggfile))

    if os.path.isfile(aggfile):
        f = file(aggfile, 'r')
        agg_output = pickle.load(f)
        f.close()
    else:
        sr, output = go_sweep_mpi(p, sweep_params, dis_class, ContactMatrix(), go_output_aggregate)
        agg_output = aggregate_over_seeds(output, p['num_runs'])
        f = file(aggfile, 'w')
        pickle.dump(agg_output, f)
        f.close()

    return agg_output


def convert_old(p, disease_type, cmatrix, cur_seed, p_string=None):
    rng = Random(cur_seed)
    disease_fname = os.path.join(p['prefix'], 'disease.hd5')

    if not os.path.isfile(disease_fname):
        print(("Problem: disease (%sseed=%d)..." % (p_string + ";" if p_string else None, cur_seed), "does not exist!"))
    else:
        disease = disease_type(p, cmatrix, rng, disease_fname, mode='a')
        disease.store_params(p)
        disease.done(True)

import shutil


def remove_incomplete_parallel(obj, fh, rank):
    print("@@@_remove_incomplete_parallel")
    remove_incomplete(obj['p'], obj['disease_type'], obj['cmatrix'], obj['seed'], p_string=obj['p_string'])


def remove_incomplete(p, disease_type, cmatrix, cur_seed, p_string=None):
    rng = Random(cur_seed)
    disease_fname = os.path.join(p['prefix'], 'disease.hd5')

    if not os.path.isfile(disease_fname):
        print(("Problem: disease (%sseed=%d)..." % (p_string + ";" if p_string else None, cur_seed), "does not exist!"))
    else:
        try:
            disease = disease_type(p, cmatrix, rng, disease_fname, mode='r')
            if disease.is_complete():
                disease.close()
            else:
                disease.close()
                if p['remove']:
                    print(("will remove", p['prefix']))
                    shutil.rmtree(p['prefix'])
        except:
            print(("Problem: disease (%sseed=%d)..." % (p_string + ";" if p_string else None, cur_seed), "exists but is empty!"))
            if p['remove']:
                print(("will remove", p['prefix']))
                shutil.rmtree(p['prefix'])

