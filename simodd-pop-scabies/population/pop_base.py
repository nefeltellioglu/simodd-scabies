"""
.. module:: pop_base
.. moduleauthor:: Nic Geard <ngeard@unimelb.edu.au>
"""

from collections import defaultdict

from population.individual import Individual


class Population(object):
    """
    The base class for a population containing objects of type 
    :class:`.individual.Individual`, or subclass.
    
    :class:`.population.Population` is the simplest population class, providing
    sufficient functionality for age-structured populations and arbitrary
    contact groups.  While it provides the basic functionality necessary for
    implementing households and demographic processes, it stops short of
    implementing these explicitly.

    A derived class, :class:`.pop_hh.PopHH`,  implements households (a type of
    contact group) and associated functionality.

    :param ind_type: the :class:`.individual.Individual` (sub)class stored by this population.
    :type ind_type: class

    .. note::

        An issue arises when an individual is added to two groups of the same
        type; their group membership is changed, but the membership list of
        their old group is NOT updated.  I need to decide how best to enforce
        this.
    
    """

    def __init__(self, ind_type=Individual):
        self.ind_type = ind_type

        # A dictionary of individuals keyed by their unique ID.
        self.I = {}

        # A dictionary of individual IDs keyed by age
        # This uses some space, but provides much more rapid access when trying
        # to select individuals by age (e.g., as mothers or partners)
        self.I_by_age = defaultdict(dict)

        # A nested dictionary containing contact group memberships.  The first
        # key specifies the type of group (a string) and the second key 
        # specifies the group's ID.  The value is a list of IDs of individuals 
        # belonging to that group.
        self.groups = {}

        # A dictionary of counters storing the next available ID for 
        # individuals and group types.
        self.next_id = {'individual': 0}

    # basic individual and group operations ###

    def add_individual(self, age=-1, sex=0, adam=False, logging=True):
        """
        Add a new individual to the population.

        :param age: The age of the newly added individual.
        :type age: int
        :param sex: The sex of the newly added individual.
        :type sex: int (0, 1)
        :param adam: ``True`` if individual is created as part of bootstrap population.
        :type adam: bool
        :param logging: ``True`` if storing logging information.
        :type logging: bool
        :returns: ind_type -- the newly added individual.
        """

        new_id = self.next_id['individual']
        self.I[new_id] = self.ind_type(new_id, age, sex, adam, logging)
        self.I_by_age[age][new_id] = self.I[new_id]
        self.next_id['individual'] += 1
        return self.I[new_id]

    def remove_individual(self, ind):
        """
        Remove an individual from the population and, by extension, from all
        of their groups.

        :param ind: The individual to be removed.
        :type ind: ind_type
        """

        for cur_group in list(ind.groups.keys()):
            self.remove_individual_from_group(cur_group, ind)
        del self.I_by_age[ind.age][ind.ID]
        self.I.pop(ind.ID)

    def init_group_type(self, group_type):
        """
        Add a new type of group to the population.

        :param group_type: The label used to denote this group type.
        :type group_type: string
        """

        assert group_type not in self.groups

        self.groups[group_type] = {}
        self.next_id[group_type] = 0

    def add_group(self, group_type, inds=None):
        """
        Add a new contact group of specified type to the population and
        return its ID.  Optionally, specify a set of individuals to be added to
        the newly created group.

        :param group_type: The type of group to add.
        :type group_type: string
        :param inds: individuals to add to the group.
        :type inds: list
        :returns: The ID of the newly created group.

        .. note::

            The group type must already have been initialised using :func:`init_group_type`.

        """

        assert group_type in self.groups

        new_id = self.next_id[group_type]
        self.groups[group_type][new_id] = []
        if inds:
            self.add_individuals_to_group(group_type, new_id, inds)
        self.next_id[group_type] += 1
        return new_id

    def remove_group(self, group_type, group_id):
        """
        Remove an existing contact group.

        :param group_type: The type of the group to be removed.
        :type group_type: string
        :param group_id: The ID of the group to be removed.
        :type group_id: int
        """

        assert group_id in self.groups[group_type]

        for cur_ind in self.groups[group_type][group_id]:
            self.I[cur_ind.ID].groups.pop(group_type)
        self.groups[group_type].pop(group_id)

    def add_individuals_to_group(self, group_type, group_id, inds):
        """
        Add individuals to an existing group.

        :param group_type: The type of the group.
        :type group_type: string
        :param group_id: The ID of the group.
        :type group_id: int
        :param inds: list of individuals to add to the group.
        :type inds: list

        .. note::
            
            No checking is done to ensure that individuals are not already a
            member of a group of this type.

        """

        assert group_id in self.groups[group_type]

        self.groups[group_type][group_id].extend(inds)
        for cur_ind in inds:
            cur_ind.groups.__setitem__(group_type, group_id)

    def remove_individual_from_group(self, group_type, ind):
        """
        Remove an individual from a group to which they belong.  Additionally,
        remove the group if it now contains no members.

        :param group_type: The type of the group the individual is to leave.
        :type group_type: str
        :param ind: The individual to be removed.
        :type ind: ind_type
        :returns: The ID of the group the individual has left.
        """

        assert group_type in ind.groups

        group_id = ind.groups[group_type]
        ind.groups.pop(group_type)
        self.groups[group_type][group_id].remove(ind)
        # remove group if now empty
        if len(self.groups[group_type][group_id]) <= 0:
            self.groups[group_type].pop(group_id)
        return group_id

    # update age of population ##########

    def age_population(self, period):
        """
        Age each individual in population by duration specified (in days).

        :param period: The number of days to age the population.
        :type period: int
        """

        birthdays = defaultdict(list)

        for ind in self.I.values():
            ind.age_days += period
            if ind.age_days >= 364:  # i has a birthday
                del self.I_by_age[ind.age][ind.ID]
                ind.age += 1
                ind.age_days %= 364
                birthdays[ind.age - 1].append(ind.ID)
                self.I_by_age[ind.age][ind.ID] = self.I[ind.ID]

            if ind.age < 0:
                print(ind.age, ind.age_days)
                print("wtf?!")
                exit()

        return birthdays

    # access individuals and groups by age, size, etc. ##########

    def ind_ids_by_age(self, min_age, max_age):
        """
        Return a list of individuals in the specified age range.

        :param min_age: The minimum age to include.
        :type min_age: int
        :param max_age: The maximum age to include.
        :type max_age: int
        :returns: a list of IDs of individuals in age range.
        """

        if min_age == max_age:
            return list(self.I_by_age[min_age].keys())
        else:
            inds = []
            for cur_age in range(min_age, max_age + 1):
                inds += list(self.I_by_age[cur_age].keys())
            return inds

    def individuals_by_age(self, min_age, max_age=None):
        """
        Return a list of individuals in the specified age range.

        :param min_age: The minimum age to include.
        :type min_age: int
        :param max_age: The maximum age to include.
        :type max_age: int
        :returns: a list of individuals in age range.
        """

        if (max_age is None) or (min_age == max_age):
            return list(self.I_by_age[min_age].values())
        else:
            inds = []
            for cur_age in range(min_age, max_age + 1):
                inds += list(self.I_by_age[cur_age].values())
            return inds

    def individuals_by_group_size(self, group_type, size):
        """
        Return a list of individuals in groups of specified type and size.

        :param group_type: The type of groups to aggregate.
        :type group_type: string
        :param size: The size of groups to aggregate.
        :type size: int
        :returns: A list of individuals in appropriately sized groups.

        """

        return [x for y in [grp for grp in list(self.groups[group_type].values())
                            if len(grp) == size] for x in y]

    def all_individuals_by_group_size(self, group_type, max_size):
        """
        Return a dictionary of individuals in groups of specified
        type, keyed by size.  All groups of max_size or larger are binned.

        :param group_type: The type of groups to aggregate.
        :type group_type: string
        :param max_size: the maximum possible group size.
        :type max_size: int
        :returns: A list of individuals in appropriately sized groups.
        """

        by_size = dict([(x, []) for x in range(1, max_size + 1)])
        for grp in self.groups[group_type].values():
            by_size[min(max_size, len(grp))].extend(grp)
        return by_size

    def individuals_by_min_group_size(self, group_type, size):
        """
        Return a list of individuals in groups of specified type and
        minimum size.

        :param group_type: The type of groups to aggregate.
        :type group_type: string
        :param size: The size of groups to aggregate.
        :type size: int
        :returns: A list of individuals in appropriately sized groups.

        """

        return [x for y in [grp for grp in list(self.groups[group_type].values())
                            if len(grp) >= size] for x in y]

    def all_groups_by_size(self, group_type, max_size):
        """
        Return a dictionary mapping size to a list of groups of that size
        """

        by_size = defaultdict(list)
        for grp in self.groups[group_type].values():
            by_size[min(max_size, len(grp))].append(grp)
        return by_size

    def groups_by_size(self, group_type, size):
        """
        Return a list of IDs of groups that are of the specified type and size.

        :param group_type: The type of groups to evaluate.
        :type group_type: string
        :param size: The size of groups to evaluate.
        :type size: int
        :returns: A list of IDs of groups of appropriate size.
        """

        return [x for x in list(self.groups[group_type].keys())
                if len(self.groups[group_type][x]) == size]

    def groups_by_min_size(self, group_type, size):
        """
        Return a list of IDs of groups that are of the specified type and 
        *at least* the specified size.

        :param group_type: The type of groups to evaluate.
        :type group_type: string
        :param size: The minimum size of groups to evaluate.
        :type size: int
        :returns: A list of IDs of groups of at least size.
        """

        return [x for x in list(self.groups[group_type].keys())
                if len(self.groups[group_type][x]) >= size]
