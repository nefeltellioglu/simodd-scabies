import pickle
import csv


def save_pop(pop, file_path):
    """
    Save a population to a pickle file.

    :param pop: The population
    :type pop: PopHH
    :param file_path: The path to save the population to.
    :type file_path: string
    """

    _flatten(pop)
    pop_file = open(file_path, 'wb')
    pickle.dump(pop, pop_file, pickle.HIGHEST_PROTOCOL)
    pop_file.close()
    _lift(pop)


def load_pop(file_path):
    """
    Load a population from a pickle file.

    :param file_path: The path to load the population from.
    :type file_path: string
    :returns: The loaded population
    """

    pop_file = open(file_path, 'rb')
    pop = pickle.load(pop_file)
    _lift(pop)
    return pop


def write_csv(pop, file_path, hh_id_map):
    """
    Write a CSV file containing population info (for RB spatial rota model)

    :param pop: The population
    :type pop: PopHH
    :param file_path: The path to save the population to
    :type file_path: string
    :param hh_id_map: map of household IDs
    :type hh_id_map: dictionary mapping simultion IDs to external IDs
    :return:
    """

    cur_id = 1
    id_map = {}

    with open(file_path, 'w') as csvfile:
        indwriter = csv.writer(csvfile)
        indwriter.writerow(['id', 'age', 'sex', 'household_id'])
        for ind in pop.I.values():
            id_map[ind.ID] = cur_id
            hh_id = hh_id_map[ind.groups['household']]
            indwriter.writerow([cur_id, ind.age, ind.sex, hh_id])
            cur_id += 1


def write_hh_csv(pop, file_path, hh_loc_path, rng):
    """
    Write a CSV file containing household info (for RB spatial rota model)

    :param pop: The population
    :type pop: PopHH
    :param file_path: The path to save the population to
    :type file_path: string
    :param hh_loc_path: Path to a file containing household locations
    :type hh_loc_path: string
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`
    :return: a map from simulation household IDs to sequential IDs
    """

    cur_hh_id = 1
    hh_id_map = {}
    hh_locs = []

    # read in household location file (generated in R)
    with open(hh_loc_path) as locfile:
        locreader = csv.reader(locfile)
        next(locreader)
        for hh in locreader:
            hh_locs.append(tuple([float(hh[1]), float(hh[2])]))

    # sample a sufficient number of household locations
    cur_locs = [hh_locs[x] for x in rng.sample(range(len(hh_locs)),
                                               len(pop.groups['household']))]

    # write households
    with open(file_path, 'w') as csvfile:
        hhwriter = csv.writer(csvfile)
        hhwriter.writerow(['id', 'latitude', 'longitude'])
        for hh in pop.groups['household'].keys():
            hh_id_map[hh] = cur_hh_id
            hhwriter.writerow([cur_hh_id, cur_locs[cur_hh_id - 1][0], cur_locs[cur_hh_id - 1][1]])
            cur_hh_id += 1

    return hh_id_map


def _flatten(pop):
    """
    For pickling... flatten recursive references.
    """

    for ind in list(pop.I.values()) + list(pop.graveyard.values()):
        if ind.partner:
            ind.partner = ind.partner.ID
        if ind.parents:
            ind.parents = [x.ID for x in ind.parents]
        if ind.children:
            ind.children = [x.ID for x in ind.children]
        if ind.deps:
            ind.deps = [x.ID for x in ind.deps]


def _lift(pop):
    """
    For unpickling... replace recursive references.
    """

    for ind in list(pop.I.values()) + list(pop.graveyard.values()):
        if ind.partner:
            if ind.partner in pop.I:
                ind.partner = pop.I[ind.partner]
            else:
                ind.partner = pop.graveyard[ind.partner]
        ind.parents = _lift_list(pop, ind.parents)
        ind.children = _lift_list(pop, ind.children)
        ind.deps = _lift_list(pop, ind.deps)


def _lift_list(pop, data):
    """
    For unpickling... replace recursive references.
    """

    lifted_data = []
    for i in data:
        if i in pop.I:
            lifted_data.append(pop.I[i])
        elif i in pop.graveyard:
            lifted_data.append(pop.graveyard[i])
    return lifted_data
