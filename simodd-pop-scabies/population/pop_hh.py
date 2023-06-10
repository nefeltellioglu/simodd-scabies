"""
.. module:: pop_hh
.. moduleauthor:: Nic Geard <ngeard@unimelb.edu.au>
"""

from population.pop_base import Population
from population.individual import Individual
from population.household import Household
from population.pop_info import hh_size


class PopHH(Population):
    """
    A population class in which *households* are a fundamental unit of organisation.  

    :param ind_type: the :class:`.individual.Individual` (sub)class stored by this population.
    :type ind_type: class
    :param logging: whether or not to write ind/hh logs
    :type logging: bool


    """

    def __init__(self, ind_type=Individual, logging=True):
        super(PopHH, self).__init__(ind_type)

        super(PopHH, self).init_group_type('household')

        self.logging = logging
        if self.logging:
            self.households = {}
        self.graveyard = {}
        self.preg_current = {}  # map of mother IDs to due time

    def birth(self, t, mother, sex=0):
        """
        Add a newborn child with specified parents to population.

        By default, the new individual's household is set to that of the first
        parent.

        :param t: The current time step.
        :type t: int
        :param mother: The mother.
        :type mother: :class:`.individual.Individual`
        :param sex: The sex of the new individual
        :type sex: int
        :returns: The new individual.

        """

        # create the new individual
        new_ind = self.add_individual(0, sex, logging=self.logging)
        new_ind.birth_order = len(mother.children) + 1

        # assign parent and dependency relationships
        mother.set_prev_birth()
        father = mother.partner
        parents = [x for x in (mother, father) if x is not None]
        new_ind.parents = parents
        for cur_parent in parents:
            cur_parent.children.append(new_ind)
            cur_parent.deps.append(new_ind)

        hh_id = parents[0].groups['household']
        self.add_individuals_to_group('household', hh_id, [new_ind])

        if self.logging:
            new_ind.household = self.households[hh_id]
            # record own birth
            new_ind.add_log(t, 'b', "%d born to %s" % (
                new_ind.ID, [x.ID for x in (mother, father) if x is not None]))
            # record children's birth
            for cur_parent in parents:
                cur_parent.add_log(
                    t, 'c', "Gave birth to %d" % new_ind.ID, new_ind.ID)
            # add the new individual to its parent's household
            self.households[hh_id].add_log(
                t, 'cb', "Birth of child", len(self.groups['household'][hh_id]))

        return new_ind

    def death(self, t, ind):
        """
        Remove individual from population, and return a list of orphaned children.

        :param t: The current time step.
        :type t: int
        :param ind: The dead individual.
        :type ind: ind_type
        :returns: a list of orphans
        """

        # identify individuals who will be orphaned by this death
        orphans = ind.deps if not ind.partner else []

        if self.logging:
            ind.add_log(t, 'd', "Died at age %d" % ind.age)
            if ind.partner:
                ind.partner.add_log(t, 'md', "Partner %d died" % ind.ID, ind.ID)
            for cur_dep in ind.deps:
                cur_dep.add_log(t, 'gd', "Parent %d died" % ind.ID, ind.ID)
                # self.graveyard[ind.ID] = ind
        
        # remove as partner
        if ind.partner:
            ind.partner.partner = None

        # remove the dead individual's guardian(s)
        self._remove_from_guardians(
            t, ind, 'cd', "Lost dependent %d (death)" % ind.ID)
        # remove dead individual from household and population 
        self._remove_individual_from_hh(t, ind, 'd', "Member died")
        self.remove_individual(ind)

        return orphans

    def form_couple(self, t, ind, partner):
        """
        Attempt to form a new couple household.  A new couple household is 
        only formed if a suitable partner can be found.

        Return a tuple containing the individual whose household the couple
        now live in (or None if it is a new household), and the household.

        :param t: the current time step.
        :type t: int
        :param ind: The individual forming a couple.
        :type ind: ind_type
        :param partner: The new partner.
        :type partner: ind_type
        """

        assert not ind.partner and not partner.partner
        assert ind.groups['household'] != partner.groups['household']

        # form couple
        ind.partner = partner
        partner.partner = ind

        if self.logging:
            ind.add_log(t, 'm', "Marriage to %d" % partner.ID, partner.ID)
            partner.add_log(t, 'm', "Marriage to %d" % ind.ID, ind.ID)

        # TODO: rather than call sep functions, use this logic to set
        # a target_hh ID variable; if function called with None, then
        # create a new household, otherwise use this (and don't add
        # existing individuals to moving list.
        if ind.with_parents:
            if partner.with_parents:
                # both individuals live at home, create a new hh
                return None, self._form_couple_hh(t, [ind, partner])
            else:  # partner has own hh, move into it
                return partner, self._merge_hh(
                    t, (ind, partner), partner.groups['household'])
        else:  # ind has own hh
            if partner.with_parents or hh_size(self, ind) > \
                    hh_size(self, partner):
                # partner lives at home, or in a smaller household than ind
                return ind, self._merge_hh(
                    t, (partner, ind), ind.groups['household'])
            else:
                # ind lives in a smaller household than partner
                return partner, self._merge_hh(
                    t, (ind, partner), partner.groups['household'])

    def separate_couple(self, t, ind):
        """
        Separate the couple involving the specified individual, moving their
        partner into a new, single-person household.  Children remain in the 
        original household.

        :param t: The current time step.
        :type t: int
        :param ind: The individual to separate.
        :type ind: ind_type
        """

        assert ind.partner

        if ind.sex == 0:
            ind_m = ind
            ind_f = ind_m.partner
        else:
            ind_f = ind
            ind_m = ind_f.partner

        ind_m.divorced = True
        ind_m.partner = None
        ind_f.divorced = True
        ind_f.partner = None
        ind_m.deps = []  # ind_f keeps kids!

        # update logs of parents and children
        if self.logging:
            ind_f.add_log(
                t, 's', "Splitting from %d, staying put" % ind_m.ID, ind_m.ID)
            ind_m.add_log(
                t, 's', "Splitting from %d, moving out" % ind_f.ID, ind_f.ID)
            for cur_dep in ind_f.deps:
                cur_dep.add_log(
                    t, 'gs', "Parents divorcing, staying with %d" % ind_f.ID)

        self._form_single_hh(t, ind_m)

    def leave_home(self, t, ind):
        """
        Move an individual out of their home, removing them as a dependent of
        their household head(s).

        :param t: the current time step.
        :type t: int
        :param ind: the individual to leave home.
        :type ind: ind_type
        """
        assert not ind.partner, "Time %s ID %s"%(t, ind.ID)
        # assert self.hh_size(ind) > 1

        ind.with_parents = False
        self._form_single_hh(t, ind)

    def allocate_orphan(self, t, ind, tgt_hh):
        """
        Move children who have been orphaned by the death of their final
        remaining parent (in the same household) to new, randomly chosen
        households.

        :param t: The current time step.
        :type t: int
        :param ind: The orphaned individual.
        :type ind: ind_type
        :param tgt_hh: The new household to move the orphan to.
        :type tgt_hh: int
        """

        assert not ind.partner
        assert not ind.children

        # remove them from their current household
        self._remove_individual_from_hh(t, ind, 'c-', "Lost relocated child")

        # appoint eldest person in new household and partner as guardians
        g_ind = sorted(self.groups['household'][tgt_hh],
                       key=lambda x: x.age)[-1]
        g_ind.deps.append(ind)
        if g_ind.partner:
            g_ind.partner.deps.append(ind)

        # add individual to new household (must happen last!!)
        self.add_individuals_to_group('household', tgt_hh, [ind])

        if self.logging:
            g_ind.add_log(t, 'c+', "Gained dependent %d (orphan)" % ind.ID)
            if g_ind.partner:
                g_ind.partner.add_log(t, 'c+',
                                      "Gained dependent %d (orphan)" % ind.ID)
            self.households[tgt_hh].add_log(
                t, 'c+', "Gained relocated child",
                len(self.groups['household'][tgt_hh]))
            ind.add_log(t, 'r', "Relocated - with %d as guardian" % g_ind.ID)

    def _form_single_hh(self, t, ind):
        """
        Move specified individual, with no partner, into their own household.

        :param t: the current time step.
        :type t: int
        :param ind: the first partner in the couple.
        :type ind: ind_type
        :returns: the ID of the new household.
        """

        assert not ind.partner
        #        assert len(self.housemates(i_id)) > 0

        # remove any guardians
        self._remove_from_guardians(
            t, ind, 'c-', "Lost dependent %d (leaving)" % ind.ID, ind.ID)

        # remove from old household 
        self._remove_individual_from_hh(t, ind, 'l', "Individual left")

        # add to newly created household
        new_hh = self.add_group('household', [ind])

        ind.with_parents = False

        if self.logging:
            self.households[new_hh] = Household(t)
            ind.add_log(t, 'l1', "Leaving household (single)")
            self.households[new_hh].add_log(t, 'f1', "Household formed", 1)

        return new_hh

    def _merge_hh(self, t, inds, hh_id):
        """
        Move p_id and any of their dependents into hh_id.

        :param t: The current time step.
        :type t: int
        :param inds: The new couple.
        :type inds: tuple
        :param hh_id: The household they will reside in.
        :param hh_id: int
        """

        new_inds = [inds[0]]
        combined_deps = inds[0].deps + inds[1].deps
        inds[0].deps = combined_deps[:]
        inds[1].deps = combined_deps[:]
        new_inds.extend(inds[1].deps)

        inds[0].with_parents = False
        inds[1].with_parents = False

        # remove any guardians of moving partner
        self._remove_from_guardians(
            t, inds[0], 'c-',
            "Lost dependent %d (leaving)" % inds[0].ID, inds[0].ID)

        for ind in new_inds:
            self._remove_individual_from_hh(t, ind, 'l', "Individual left")

        self.add_individuals_to_group('household', hh_id, new_inds)

        if self.logging:
            inds[0].add_log(
                t, 'l2',
                "Leaving household - couple with %d" % inds[1].ID, inds[1].ID)
            self.households[hh_id].add_log(
                t, 'm',
                "Household merged (%d individuals)" % len(new_inds),
                len(self.groups['household'][hh_id]))

        return hh_id

    def _form_couple_hh(self, t, inds):
        """
        Make specified single individuals into a couple and move them into a 
        new household.  Any dependents of either individual accompany them to 
        the new household.

        :param t: the current time step.
        :type t: int
        :param inds: the individuals to move into a couple household.
        :type inds: list
        :returns: the ID of the new household.
        """

        ind_a = inds[0]
        ind_b = inds[1]

        inds[0].with_parents = False
        inds[1].with_parents = False

        # move dependents along with guardians
        new_inds = list(inds)
        combined_deps = ind_a.deps + ind_b.deps
        ind_a.deps = combined_deps[:]
        ind_b.deps = combined_deps[:]
        new_inds.extend(combined_deps)

        # remove any guardians of new couple
        self._remove_from_guardians(
            t, inds[0], 'c-',
            "Lost dependent %d (leaving)" % inds[0].ID, inds[0].ID)
        self._remove_from_guardians(
            t, inds[1], 'c-',
            "Lost dependent %d (leaving)" % inds[1].ID, inds[1].ID)

        # remove individuals from prior households and create new household
        for ind in new_inds:
            self._remove_individual_from_hh(t, ind, 'l', "Individual left")
        hh_id = self.add_group('household', new_inds)

        if self.logging:
            self.households[hh_id] = Household(t)
            ind_a.add_log(
                t, 'l2', "Leaving household - couple with %d" % ind_b.ID)
            ind_b.add_log(
                t, 'l2', "Leaving household - couple with %d" % ind_a.ID)
            self.households[hh_id].add_log(
                t, 'f2',
                "Household formed (%d individuals)" % len(new_inds),
                len(new_inds))

        return hh_id

    def _remove_from_guardians(self, t, ind, log_code, log_entry, other=None):
        """
        Remove i_id as a dependent upon any of the other individuals in their
        household.  Appends a custom entry to the former guardian's log.

        :param ind: The individual leaving their guardians.
        :type ind: ind_type
        :log_code: The log code classification why the individual is leaving.
        :type log_code: string
        :log_msg: The log message describing why the individual is leaving.
        :type log_entry: string
        :returns: `True` if ind was dependent upon anyone.
        """

        was_dep = False
        for cur_hmate in self.groups['household'][ind.groups['household']]:
            if ind in cur_hmate.deps:
                cur_hmate.deps.remove(ind)
                if self.logging:
                    cur_hmate.add_log(t, log_code, log_entry, other)
                was_dep = True
        return was_dep

    def _remove_individual_from_hh(self, t, ind, log_code, log_msg):
        """
        Remove an individual from a household.
        """

        # TODO: should this handle dependent/guardian links?

        old_hh = self.remove_individual_from_group('household', ind)
        new_size = 0 if old_hh not in self.groups['household'] \
            else len(self.groups['household'][old_hh])

        if self.logging:
            self.households[old_hh].add_log(t, log_code, log_msg, new_size)



