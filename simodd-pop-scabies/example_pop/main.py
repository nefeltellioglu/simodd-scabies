"""
A stub for testing/evaluating dynamic demography simulations.
"""
import sys

#sim_path = '/home/unimelb.edu.au/ngeard/PycharmProjects/simodd-pop'
#sys.path.append(sim_path)

from random import Random

from population.simulation import Simulation
from observers.obs_pop import PopulationObserver
from pop_explore.output_pop import *
from example_pop.params import p


def run_single(p):
    if p['random_seed']:
        p['seed'] = Random().randint(0, 99999999)

    print("Creating population...")
    #pop_fname = os.path.join(p['prefix'], 'population.hd5')
    sim = Simulation(p)
    sim.add_observers(PopulationObserver(sim.h5file))
    iterations = p['years'] * p['t_per_year']

    print("Running simulation...")
    print("iter\tyears\tdays\tpeople\thouses\tbirths\tdeaths\tmigrants")
    for i in range(iterations):
        t = i * 364 // p['t_per_year']
        b, d, im, bd = sim.update_all_demo(i)
        sim.update_observers(t, pop=sim.P)
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            i, t // 364, t % 364,
            len(sim.P.I), len(sim.P.groups['household']),
            len(b), len(d), len(im)))

    print("Simulation complete!")
    sim.done(complete=True)
    return sim


def go_single(p):
    if p['random_seed']:
        p['seed'] = Random().randint(0, 99999999)
    sim = run_single(p)
    return sim


if __name__ == '__main__':
    cur_params = p
    cur_params['ext'] = 'svg'
    cur_params['preg'] = False
    cur_sim = run_single(cur_params)


