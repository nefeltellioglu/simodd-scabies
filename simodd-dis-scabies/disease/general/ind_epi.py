"""
Base individual class for epidemic simulations, adding state and counters.
"""
from population.individual import Individual


class IndEpi(Individual):
    __slots__ = 'prev_state', 'state', 'next_state', \
                'age_at_infection', 'vaccines_received', 'vaccines_effective', 'vaccine_times', \
                'infections', 'counters', 'source', 'hh_frac', \
                'hh_source', 'durations', 'birth_state'

    def __init__(self, new_id, age=0, sex=0, bootstrap=False, logging=True):
        super(IndEpi, self).__init__(new_id, age, sex, bootstrap, logging)
        self.prev_state = None
        self.state = None
        self.next_state = None

        self.birth_state = None

        # store age at infection to ensure infected-by-age counts are 
        # handled correctly when an individual ages WHILE infected.
        self.age_at_infection = None
        self.vaccines_received = []  # a list of vaccine types received
        self.vaccines_effective = []  # a corresponding list of Boolean values indicating whether vaccine was effective
        self.vaccine_times = []  # a list of vaccination times
        self.infections = []  # a list of infection times
        self.counters = {}

        self.source = -1  # ID of the individual responsible for most recent infection
        self.hh_frac = 0  # fraction of most recent infection attributable to hh sources
        self.hh_source = False
        self.durations = {}  # dictionary of durations in each state (key)

