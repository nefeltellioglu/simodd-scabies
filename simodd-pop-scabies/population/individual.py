"""
.. module:: individual
.. moduleauthor:: Nic Geard <ngeard@unimelb.edu.au>
"""


class Individual(object):
    """
    The base class for an individual, storing information on demographic
    attributes, population inks (family, household, etc.) .

    This class will typically be inherited to add disease, etc., information.

    :param new_id: a unique identifier.
    :type new_id: int
    :param age: the individual's age in years.
    :type age: int
    :param sex: the individual's sex (0 == male; 1 == female)
    :type sex: int
    :param adam: `True` if individual is part of the bootstrapped population.
    :type adam: bool
    :param logging: whether or not to write log
    :type logging: bool

    """

    __slots__ = 'ID', 'sex', 'age', 'age_days', 'prev_birth_age', \
                'next_birth_age', 'bootstrap', 'log', 'birth_order', 'groups', \
                'partner', 'divorced', 'parents', \
                'children', 'deps', 'with_parents'

    def __init__(self, new_id, age=0, sex=0, bootstrap=False, logging=True):
        #: unique identifier for the individual.
        self.ID = new_id
        #: sex of the individual (0 == female).
        self.sex = sex
        #: age of the individual in years.
        self.age = age
        #: fractional age of the individual in days (0--364).
        self.age_days = 0
        #: age at which individual last gave birth.
        self.prev_birth_age = 0
        #: age at which individual is next eligible to give birth.
        self.next_birth_age = 0
        #: flag indicating individual is part of a bootstrap population.
        self.bootstrap = bootstrap
        #: index of birth order (to mother).
        self.birth_order = 0

        #: dictionary of group_type:group_ID pairs
        self.groups = {}

        #: the partner ID of this individual (if any).
        self.partner = None
        #: `True` if this individual has been divorced.
        self.divorced = False
        #: the parent IDs of this individual (if known).
        self.parents = []
        #: the children IDs of this individual (if any).
        self.children = []
        #: the dependant IDs of this individual (if any).
        self.deps = []

        # NB: need this if want to store (& compare) hh_type over time
        # self.hh_type = []

        #: `True` if individual is living with one or more parents
        self.with_parents = True

        if logging:
            #: (optionally) a link to a Household object (used to store hh data and logs)
            self.household = None

            #: (optionally) a list of significant events in this individual's life.
            self.log = []

    def set_prev_birth(self):
        """
        Store the age at which birth has occurred.
        """
        self.prev_birth_age = self.age * 364 + self.age_days
        # NB: calculating minimum gap here is more efficient than in can_birth
        # as only created a single time (per birth)
        self.next_birth_age = self.prev_birth_age + 300

    def can_birth(self):
        """
        Returns `True` if at least nine months have passed since previous birth.

        :returns: `True` if individual can give birth.
        """
        return self.age * 364 + self.age_days - self.next_birth_age > 0

    def add_log(self, time, code, msg, other_id=None):
        """
        Add a log entry.

        :param time: current time step
        :type time: int
        :param code: a code for filtering events.
        :type code: string
        :param msg: a text message describing this event.
        :type msg: string
        :param other_id: ID of other individual (if applicable)
        :type other_id: int

        Possible log codes for individuals:

        - b: own birth
        - cb: birth of child
        - cd: death of child (dependent)
        - c+: gain of child (relocated)
        - c-: child leaving home
        - d: own death
        - gd: death of parent (guardian)
        - gs: parents (guardians) separating
        - l1: leaving home (single)
        - l2: leaving home (couple)
        - m: marriage (couple formation)
        - s: couple separation
        - r: relocation (as orphan)
        """

        self.log.append({
            'age': self.age,
            'code': code,
            'msg': msg,
            'other': other_id,
            'time': time
        })
