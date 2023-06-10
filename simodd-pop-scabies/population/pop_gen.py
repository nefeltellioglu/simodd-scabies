# Population generation functions #######################################
import itertools

from population.pop_hh import Household
from population.utils import sample_table, split_age_probs


def duplicate_household(pop, t, hh_id):
    """
    Create a duplicate of household hh_id, consisting of
    new individuals with same household age composition.

    Return id of newly created household.

    :param pop: The population
    :type pop: PopHH
    :param t: The current time step.
    :type t: int
    :param hh_id: The household to duplicate
    :type hh_id: int
    :returns: the id of the newly created household.
    """

    new_hh = []
    for ind in pop.groups['household'][hh_id]:
        new_ind = pop.add_individual(ind.age, ind.sex,
                                     logging=pop.logging)
        if pop.logging:
            new_ind.add_log(t, 'i', "Individual (immigrated)")
        new_hh.append(new_ind)
    new_hh_id = pop.add_group('household', new_hh)
    init_household(pop, new_hh_id)
    if pop.logging:
        pop.households[new_hh_id] = Household(t)
        pop.households[new_hh_id].add_log(
            t, 'i', "Household (immigrated)",
            len(pop.groups['household'][new_hh_id]))

    return new_hh_id


def gen_hh_size_structured_pop(pop, pop_size, hh_probs, rng):
    """
    Generate a population of individuals with household size structure.

    :param pop: The population
    :type pop: PopHH
    :param pop_size: The size of the population to be generated.
    :type pop_size: int
    :param hh_probs: A table of household size probabilities.
    :type hh_probs: list
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`
    """

    i = 0
    while i < pop_size:
        size = int(sample_table(hh_probs, rng)[0])
        cur_hh = []
        for _ in itertools.repeat(None, size):
            cur_ind = pop.add_individual(logging=pop.logging)
            cur_hh.append(cur_ind)
            i += 1
        pop.add_group('household', cur_hh)


def gen_hh_age_structured_pop(pop, pop_size, hh_probs, age_probs_i,
                              cutoffs, rng):
    """
    Generate a population of individuals with age structure and household 
    composition.

    Household composition here is approximated by the number of 
    individuals who are:

    - pre-school age (0--4)
    - school age (5--17)
    - adult (18+)

    This is a bit ad hoc, but serves current purposes, and populations are
    intended to be 'burnt in' to create more natural structures.

    :param pop: The population
    :type pop: PopHH
    :param pop_size: The size of the population to be generated.
    :type pop_size: int
    :param hh_probs: A table of household size probabilities.
    :type hh_probs: list
    :param age_probs_i: A table of age probabilities.
    :type age_probs_i: list
    :param cutoffs: A tuple of values that specify ages cutoffs between pre-school, school and adults.
    :type cutoffs: tuple
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`
    """

    age_probs = [[x, int(y[0])] for x, y in age_probs_i]
    split_probs = split_age_probs(age_probs, cutoffs)
    split_probs.reverse()

    i = 0
    while i < pop_size:
        # get list of [adults, school age, preschool age]
        hh_type = [int(x) for x in sample_table(hh_probs, rng)]
        #            hh_type = [int(x) for x in sample_uniform(hh_probs, rng)]
        cur_hh = []
        for cur_hh_type, cur_prob in zip(hh_type, split_probs):
            for _ in itertools.repeat(None, cur_hh_type):
                sex = 0 if rng.random() < 0.5 else 1
                cur_age = sample_table(cur_prob, rng)
                cur_ind = pop.add_individual(
                    cur_age, sex, adam=True, logging=pop.logging)
                if pop.logging:
                    cur_ind.add_log(0, 'f', "Individual (bootstrap)")
                cur_hh.append(cur_ind)
                i += 1
        hh_id = pop.add_group('household', cur_hh)

        if pop.logging:
            pop.households[hh_id] = Household(0, adam=True)
            for cur_ind in cur_hh:
                cur_ind.household = pop.households[hh_id]
            pop.households[hh_id].add_log(
                0, 'f', "Household (bootstrap)",
                len(pop.groups['household'][hh_id]))


def gen_single_hh_pop(pop, hh_size, rng):
    """
    Generate a dummy population consisting of a single household.

    :param pop: The population
    :type pop: PopHH
    :param hh_size: The size of the single household.
    :type hh_size: int
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`
    """
    cur_hh = []
    for _ in range(hh_size):
        sex = 0 if rng.random() < 0.5 else 1
        cur_ind = pop.add_individual(20, sex,
                                     adam=True, logging=pop.logging)
        if pop.logging:
            cur_ind.add_log(0, 'f', "Individual (bootstrap)")
        cur_hh.append(cur_ind)
    hh_id = pop.add_group('household', cur_hh)

    if pop.logging:
        pop.households[hh_id] = Household(0, adam=True)
        pop.households[hh_id].add_log(0, 'f', "Household (bootstrap)",
                                      len(pop.groups['household'][hh_id]))


def construct_hh_age_pop(pop, age_list):
    """
    Generate a population from a list of households (lists of ages).
    Currently, this is just for testing purposes, so all individuals
    are created female.

    :param pop: The population
    :type pop: PopHH
    :param age_list: a list of lists (households of ages)
    """
    for cur_ages in age_list:
        cur_hh = []
        for cur_age in cur_ages:
            cur_ind = pop.add_individual(cur_age, 0, adam=True, logging=pop.logging)
            cur_hh.append(cur_ind)
        hh_id = pop.add_group('household', cur_hh)

        if pop.logging:
            pop.households[hh_id] = Household(0, adam=True)
            for cur_ind in cur_hh:
                cur_ind.household = pop.households[hh_id]


def allocate_couples(pop):
    """
    For testing/bootstrapping: given a household-structured population,
    for all households containing two or more people who are 18 or older,
    form the two oldest members of the household into a couple and 
    make remaining members dependents of the household head(s).

    There is a rather ugly hack here at the moment to prevent a couple
    being created with a dependent who is the same age.  Initial hh 
    allocation possibly needs to be a touch more elegant.

    :param pop: The population
    :type pop: PopHH
    """

    for cur_hh in list(pop.groups['household'].keys()):
        init_household(pop, cur_hh)


def init_household(pop, cur_hh):
    """
    Bootstrap initial household relationships; making two eldest individuals
    partners (if > 17 years).

    Modifies ages if having two individuals of equal age is likely to cause
    problems (usually for reallocate orphans).

    Again, a fairly ugly hack at the moment.

    :param pop: The population
    :type pop: PopHH
    :param cur_hh: The household to initialise
    :type cur_hh: The identifier of the household
    """

    hh_size = len(pop.groups['household'][cur_hh])
    if hh_size == 1:
        pop.groups['household'][cur_hh][0].with_parents = False
    # print "@ household of size 1"
    #            print "#0 (%d)"%(pop.groups['household'][cur_hh][0].age), \
    #                    pop.groups['household'][cur_hh][0].with_parents
    else:
        sorted_by_age = sorted(pop.groups['household'][cur_hh],
                               key=lambda x: x.age, reverse=True)

        if hh_size > 2 and sorted_by_age[1].age == sorted_by_age[2].age:
            del pop.I_by_age[sorted_by_age[2].age][sorted_by_age[2].ID]
            sorted_by_age[2].age -= 1
            if sorted_by_age[2].age < 0:
                sorted_by_age[2].age += 2
            pop.I_by_age[sorted_by_age[2].age][sorted_by_age[2].ID] = \
                sorted_by_age[2]
            if pop.logging:
                sorted_by_age[2].log[0]['age'] = sorted_by_age[2].age
        if sorted_by_age[0].age == sorted_by_age[1].age:
            del pop.I_by_age[sorted_by_age[0].age][sorted_by_age[0].ID]
            sorted_by_age[0].age += 1
            pop.I_by_age[sorted_by_age[0].age][sorted_by_age[0].ID] = \
                sorted_by_age[0]
            if pop.logging:
                sorted_by_age[0].log[0]['age'] = sorted_by_age[0].age
        if sorted_by_age[1].age > 17:
            form_couple_no_hh(sorted_by_age[:2])
            sorted_by_age[0].deps = sorted_by_age[2:]
            sorted_by_age[1].deps = sorted_by_age[2:]
            sorted_by_age[0].with_parents = False
            sorted_by_age[1].with_parents = False
        else:
            sorted_by_age[0].deps = sorted_by_age[1:]
            sorted_by_age[0].with_parents = False

            #            print "@ household of size", len(sorted_by_age)
            #            for i, xx in enumerate(sorted_by_age):
            #                print "#%d -- %d -- (%d)"%(i, xx.ID, xx.age), xx.with_parents,
            #                if xx.partner:
            #                    print xx.partner.ID
            #                else:
            #                    print ""


def form_couple_no_hh(inds):
    """
    For testing/bootstrapping: forms a couple (as above), but sets 
    ind.partner fields without modifying households.

    :param inds: The individuals to couple
    :type inds: tuple of two Individuals
    """

    assert inds[0].groups['household'] == inds[1].groups['household']

    inds[0].partner = inds[1]
    inds[1].partner = inds[0]


# initialisation functions for setting up age and group structures ###

def gen_age_structured_pop(pop, pop_size, age_probs, rng):
    """
    Generate a population of individuals with given age structure.

    :param pop: The population
    :type pop: Population
    :param pop_size: The number of individuals to generate.
    :type pop_size: int
    :param age_probs: A table mapping probabilities to age.
    :type age_probs: list
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`
    """

    for _ in itertools.repeat(None, pop_size):
        pop.add_individual(age=int(sample_table(age_probs, rng)[0]))


def allocate_groups_by_age(pop, group_type,
                           size_probs, age_lims, rng):
    """
    Allocate individuals in a given age range to groups
    with a given size distribution.

    :param pop: The population
    :type pop: Population
    :param group_type: The type of groups to create.
    :type group_type: string
    :param size_probs: A table mapping probability to group size.
    :type size_probs: list
    :param age_lims: A tuple of the min/max ages to include.
    :type age_lims: tuple
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`

    """

    assert group_type not in pop.groups

    pop.init_group_type(group_type)

    inds = pop.individuals_by_age(age_lims[0], age_lims[1])
    rng.shuffle(inds)
    while len(inds) > 0:
        size = int(sample_table(size_probs, rng)[0])
        members = [x for x in inds[:size]]
        group_id = pop.add_group(group_type, members)
        for cur_ind in members:
            cur_ind.groups.__setitem__(group_type, group_id)
        inds = inds[size:]


def allocate_groups_by_group(pop, target_type, source_type,
                             size_probs, rng):
    """
    Create 'groups of groups' by randomly aggregating lower level groups
    into higher level groups according to the given size distribution.

    For example, build neighbourhoods out of households by
    grouping households according to some distribution over number of
    households per neighbourhood.

    :param pop: The population
    :type pop: Population
    :param target_type: The type of groups to create.
    :type target_type: str
    :param source_type: The type of groups to aggregate.
    :type source_type: int
    :param size_probs: A table mapping probability to group size.
    :type size_probs: list
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`

    """

    assert source_type in pop.groups
    assert target_type not in pop.groups

    pop.init_group_type(target_type)

    ids = list(pop.groups[source_type].keys())
    rng.shuffle(ids)
    while len(ids) > 0:
        size = int(sample_table(size_probs, rng)[0])
        members = []
        group_id = pop.add_group(target_type, members)
        for source_id in ids[:size]:
            for cur_ind in pop.groups[source_type][source_id]:
                cur_ind.groups[target_type] = group_id
                members.append(cur_ind)
            pop.add_individuals_to_group(target_type, group_id, members)
        del ids[:size]
