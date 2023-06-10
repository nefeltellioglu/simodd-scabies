# Info and helper functions #############################################
import numpy
import numpy as np


def num_parents_in_hh(pop, ind):
    """
    Get the number of parents present in an individual's household.

    :param pop: The population
    :type pop: PopHH
    :param ind: The individual in question.
    :type ind: ind_type
    :returns: The number of parents.
    """
    return_value = 0
    for i in pop.hh_members(ind):
        if ind in i.deps:
            return_value += 1
    return return_value


def hh_type_stats(pop):
    """
    Aggregates data on the type of household each individual belongs to.

    Current possible household types are: single only, couple only,
    single with children, couple with children, and with parents

    Future types should include (at least, maybe) group household

    :returns: A dictionary mapping household type to a list of individuals
            belonging to a household of that type.

    :param pop: The population
    :type pop: PopHH
    """

    stats = {'single_only': [], 'couple_only': [],
             'single_kids': [], 'couple_kids': [],
             'with_parents': []}

    for i_id, ind in pop.I.items():
        cur_hh = pop.groups['household'][ind.groups['household']]
        if len(cur_hh) == 0:
            print("HH of size 0!")
        if len(cur_hh) == 1:
            stats['single_only'].append(i_id)
        else:
            if ind.with_parents:
                stats['with_parents'].append(i_id)
            else:
                if ind.partner in set(cur_hh):
                    if len(cur_hh) == 2:
                        stats['couple_only'].append(i_id)
                    else:
                        stats['couple_kids'].append(i_id)
                else:
                    stats['single_kids'].append(i_id)

    return stats


def get_hh_type(pop, ind):
    """
    Get the type of household an individual belongs to.

    :param pop: The population
    :type pop: PopHH
    :param ind: The individual in question.
    :type ind: ind_type
    :returns: The household type.
    """

    cur_hh = pop.hh_members(ind)
    if len(cur_hh) == 0:
        print("HH of size 0!")
    if len(cur_hh) == 1:
        return 'single_only'
    else:
        num_parents = pop.num_parents_in_hh(ind)
        if num_parents == 1:
            return 'single_kids'
        elif num_parents == 2:
            return 'couple_kids'
        else:
            if ind.partner in set(cur_hh):
                if len(cur_hh) == 2:
                    return 'couple_only'
                else:
                    return 'couple_kids'
            else:
                return 'single_kids'


def sum_hh_stats_group(pop):
    """
    Get a summary of the number of households of each type in population.

    :param pop: The population
    :type pop: PopHH
    :returns: A dictionary mapping household types to counts.
    """

    hh_stats = {'couple_kids': 0, 'couple_only': 0,
                'single_kids': 0, 'single_only': 0}

    for k, cur_hh in list(pop.groups['household'].items()):
        hh_type = 'with_parents'
        for cur_ind in cur_hh:
            cur_type = pop.get_hh_type(cur_ind)
            if cur_type != 'with_parents':
                hh_type = cur_type

        hh_stats[hh_type] += 1
    del hh_stats['single_only']
    fam_count = sum([v for v in list(hh_stats.values())])
    for k, v in list(hh_stats.items()):
        hh_stats[k] = float(v) / fam_count

    return hh_stats


def housemates(pop, ind):
    """
    Returns a list of ids of other individuals in i_id's household.

    :param pop: The population
    :type pop: PopHH
    :param ind: The individual to test.
    :type ind: ind_type
    :returns: A list of IDs of other individuals sharing i_id's house.
    """

    return [x for x in
            pop.groups['household'][ind.groups['household']]
            if x is not ind]


def hh_members(pop, ind):
    """
    Get a list of individuals in i_id's household (including i_id).

    :param pop: The population
    :type pop: PopHH
    :param ind: The individual to test.
    :type ind: ind_type
    :returns: The size of the individual's household.
    """

    return pop.groups['household'][ind.groups['household']]


def hh_size(pop, ind):
    """
    Get the size of an individual's household.

    :param pop: The population
    :type pop: PopHH
    :param ind: The individual to test.
    :type ind: ind_type
    :returns: The size of the individual's household.
    """

    return len(pop.groups['household'][ind.groups['household']])


def hh_age(pop, t, ind):
    """
    Get the age of an individual's household.

    :param pop: The population
    :type pop: PopHH
    :param t: The current time.
    :type t: int
    :param ind: The individual to test.
    :type ind: ind_type
    :returns: The age of the individual's household.
    """

    if pop.logging:
        return t - pop.households[ind.groups['household']].founded
    else:
        return -1


def hh_parents_siblings(pop, ind):
    """
    Get the number of parents/siblings in an individual's household.

    :param pop: The population
    :type pop: PopHH
    :param ind: The individual in question.
    :type ind: ind_type
    :returns: A tuple containing (number of parents, number of children).
    """

    num_parents = pop.num_parents_in_hh(ind)
    num_children = 0 if num_parents is 0 else \
        pop.hh_size(ind) - num_parents - 1
    return num_parents, num_children


def create_age_map(cutoffs):
    """
    Create a map from age to age bin given specified bin cutoffs.
    """

    age_map = []
    prev_cutoff = 0
    for index, cur_cutoff in enumerate(cutoffs):
        age_map.extend([index] * (cur_cutoff - prev_cutoff))
        prev_cutoff = cur_cutoff
    return age_map


def get_mean_housemates_by_age(pop, cutoffs=None, max_age=101):
    """
    Get the mean number of housemates that individuals of a particular age have.
    (Used for adjusting activity levels in community mixing matrix)

    :param pop: The population
    :type pop: PopHH
    :param cutoffs: List of age bracket cutoffs
    :type cutoffs: List of integers
    :param max_age: Upper bound on ages
    :type max_age: int
    """
    if not cutoffs:
        cutoffs = list(range(max_age))
    age_map = create_age_map(cutoffs[1:] + [max_age])
    housemate_counts = [[] for _ in range(len(cutoffs))]
    for ind in list(pop.I.values()):
        housemate_counts[age_map[ind.age]].append(
            len(pop.groups['household'][ind.groups['household']]) - 1)

    return [np.mean(x) for x in housemate_counts]


def ancestors(pop, base_ind, cur_ind=None, cur_list=None):
    """
    Return a list of all of an individual's ancestors NOT in the same household

    :param pop: The population
    :type pop: PopHH
    :param base_ind: The individual in question.
    :type base_ind: ind_type
    :param cur_ind: The current individual being considered.
    :type cur_ind: ind_type
    :param cur_list: The current list of ancestors.
    :type cur_list: list
    :returns: A list of ancestors.
    """

    if cur_list is None:
        cur_list = []

    if cur_ind is None:
        cur_ind = base_ind
        # print "base_ind = %d (%d)" % (base_ind.ID, base_ind.age)
        # print "cur_ind = %d (%d)" % (cur_ind.ID, cur_ind.age)
        # print "cur_list = %s" % ", ".join(['%d'%x.ID for x in cur_list])
        # print "parents = %s" % ", ".join(['%d'%x.ID for x in cur_ind.parents])

    for cur_anc in cur_ind.parents:
        # print "cur_anc = %d" % cur_anc.ID
        # check current ancestor is alive
        if cur_anc.ID not in pop.I:
            continue
        # print "  (is alive)"
        # check current ancestor is not in base_ind's household
        if cur_anc.groups['household'] != base_ind.groups['household']:
            # print "  (is not in base household)"
            # otherwise add current ancestor to list
            cur_list.append(cur_anc)
        # print "  (added to list... call recursively)"
        # check their ancestors
        pop.ancestors(base_ind, cur_anc, cur_list)

    return cur_list


def age_dist(I, num_bins=101, max_age=101, norm=True):
    """
    Return the age distribution of a population.

    :param num_bins: the number of bins to group the population into.
    :type num_bins: int
    :param max_age: the maximum possible age.
    :type max_age: int
    :param norm: return proportions if `True`, otherwise counts.
    :type norm: bool
    :returns: a tuple containing a list of values and a list of bin edges.
    """

    # TODO: factor out (replace with get_age_list())

    ages = [i.age for i in list(I.values())]
    return np.histogram(ages, bins=num_bins, range=(0, max_age),
                        normed=norm)


def age_dist_month(I, bin_days=None, max_age=5, norm=True):
    """
    Return the age distribution of the population (up to max_age) in bins
    of specified size.

    :param bin_days: the number of days per age bin.
    :type bin_days: int
    :param max_age: the maximum possible age.
    :type max_age: int
    :param norm: return proportions if `True`, otherwise counts.
    :type norm: bool
    :returns: a tuple containing a list of values and a list of bin edges.
    """

    if not bin_days:
        bin_days = [0, 30, 61, 89, 120, 150, 181, 211, 242, 273, 303, 334]
    age_bins = []
    for cur_age in range(max_age):
        age_bins.extend([(cur_age) * 364 + x for x in bin_days])
    age_bins.extend([364*max_age])
    # get ages in days
    ages = [i.age * 364 + i.age_days for i in list(I.values())]
    return np.histogram(ages, bins=age_bins, normed=norm)


def group_size_dist(groups, group_type, max_size=10, norm=True):
    """
    Return the size distribution of groups of specified type.

    :param group_type: the type of group to evaluate.
    :type group_type: string
    :param max_size: the maximum possible group size.
    :type max_size: int
    :param norm: return proportions if `True`, otherwise counts.
    :type norm: bool
    :returns: a tuple containing a list of values and a list of bin edges.
    """

    # TODO: factor out (replace with get_group_size_list())

    sizes = [len(hh) for hh in list(groups[group_type].values())]
    return np.histogram(sizes, bins=max_size, range=(1, max_size + 1),
                        normed=norm)


def group_size_avg(groups, group_type):
    """
    Return the average size of groups of specified type.

    :param group_type: the type of group to evaluate.
    :type group_type: string
    :returns: The mean size of groups of type `group_type`.
    """

    sizes = [len(hh) for hh in list(groups[group_type].values())]
    return np.mean(sizes)


def individuals_by_group_size_dist(groups, group_type, max_size=10,
                                   norm=True):
    """
    Return the distribution of number of individuals by group size.

    :param group_type: the type of group to evaluate.
    :type group_type: string
    :param max_size: the maximum possible group size.
    :type max_size: int
    :param norm: return proportions if `True`, otherwise counts.
    :type norm: bool
    :returns: a tuple containing a list of values and a list of bin edges.
    """

    size_dist, bins = group_size_dist(groups, group_type, max_size, norm)
    dist = [i * x for i, x in enumerate(size_dist)]
    return dist, bins